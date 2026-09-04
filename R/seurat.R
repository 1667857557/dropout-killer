#' Run DropoutKiller on a Seurat object
#'
#' In the default `modality = "rna"` path this function preserves the historical
#' behavior. With `modality = "multiome"`, RNA and ATAC share one frozen WNN
#' membership. RNA keeps the existing detector/recovery engine while ATAC uses
#' the accessibility/capture model and selective Gamma-Poisson recovery in
#' `R/atac.R`. Raw assays and fragment records are never overwritten.
#'
#' @export
dropout_killer_seurat <- function(
    object, assay = NULL, slot = "counts", reduction = "pca",
    dims = 1:20, group_by = NULL, split_by = NULL,
    new_assay = "DropoutKiller", return_result = FALSE,
    normalize = NULL, normalization_scale_factor = 1e4,
    modality = c("rna", "multiome"),
    atac_assay = "ATAC", atac_slot = "counts",
    wnn_graph = "wsnn", wnn_dims = 30L,
    rna_new_assay = NULL, atac_new_assay = "ATAC_DK",
    atac_bfdr = 0.01, atac_pi_min = 0.80,
    atac_pre_pi = 0.50, atac_score_min = 0.50,
    atac_min_observed_donors = 2L,
    atac_min_effective_donors = 5,
    atac_capture_max_iter = 25L,
    atac_capture_tol = 1e-5, atac_capture_eps = 1e-4,
    atac_max_candidates = 5e6,
    atac_phi_prior = 1, atac_phi_kappa = 10,
    atac_phi_floor = 1e-4, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required", call. = FALSE)
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object", call. = FALSE)
  modality <- match.arg(modality)
  if (is.null(assay)) {
    assay <- if (modality == "multiome" && "RNA" %in% names(object@assays)) {
      "RNA"
    } else {
      Seurat::DefaultAssay(object)
    }
  }
  meta <- object[[]]

  group <- NULL; split <- NULL
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(meta)) stop("group_by column not found", call. = FALSE)
    group <- meta[[group_by]]; names(group) <- rownames(meta)
  }
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(meta)) stop("split_by column not found", call. = FALSE)
    split <- meta[[split_by]]; names(split) <- rownames(meta)
  }

  if (modality == "rna") {
    x <- Seurat::GetAssayData(object = object, assay = assay, slot = slot)
    emb <- Seurat::Embeddings(object, reduction = reduction)
    dims <- as.integer(dims); dims <- dims[dims >= 1L & dims <= ncol(emb)]
    if (!length(dims)) stop("no requested dimensions exist in reduction", call. = FALSE)
    emb <- emb[, dims, drop = FALSE]
    if (is.null(normalize)) normalize <- identical(slot, "counts")
    if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
      stop("normalize must be TRUE or FALSE", call. = FALSE)
    if (identical(slot, "counts") && !normalize)
      warning("slot='counts' with normalize=FALSE leaves raw counts on the recovery scale", call. = FALSE)
    if (!identical(slot, "counts") && normalize)
      warning("normalizing a non-count Seurat slot; set normalize=FALSE if the selected slot is already library/log normalized", call. = FALSE)
    res <- dropout_killer(x = x, embedding = emb, group = group, split_by = split,
      normalize = normalize, normalization_scale_factor = normalization_scale_factor, ...)
    object[[new_assay]] <- Seurat::CreateAssayObject(data = res$expression)
    misc <- object@misc
    if (is.null(misc$DropoutKiller)) misc$DropoutKiller <- list()
    misc$DropoutKiller[[new_assay]] <- list(settings = res$settings,
      n_memberships = length(unique(res$membership)), n_detected = nrow(res$events),
      n_recovered = if (nrow(res$events)) sum(res$events$changed) else 0L)
    object@misc <- misc
    if (return_result) list(object = object, result = res) else object
  } else {
    if (!assay %in% names(object@assays)) stop("RNA assay not found in Seurat object", call. = FALSE)
    if (!atac_assay %in% names(object@assays)) stop("atac_assay not found in Seurat object", call. = FALSE)
    if (!wnn_graph %in% names(object@graphs))
      stop("wnn_graph not found in Seurat object@graphs; run FindMultiModalNeighbors first", call. = FALSE)

    dots <- list(...)
    take <- function(name, default) if (!is.null(dots[[name]])) dots[[name]] else default
    supplied_membership <- dots$membership; dots$membership <- NULL
    for (nm0 in c("x", "embedding", "group", "split_by", "normalize", "normalization_scale_factor"))
      dots[[nm0]] <- NULL

    graph <- object@graphs[[wnn_graph]]
    if (is.null(supplied_membership)) {
      membership_fit <- build_supercell_membership(
        graph = graph, graph_dims = wnn_dims, group = group, split_by = split,
        gamma = take("gamma", 150), k_knn = take("k_knn", 5L),
        approximate = take("approximate", "auto"), approx_n = take("approx_n", 20000L),
        seed = take("seed", 12345L))
      membership <- membership_fit$membership
    } else if (inherits(supplied_membership, "DropoutKillerMembership")) {
      membership_fit <- supplied_membership; membership <- supplied_membership$membership
    } else {
      membership_fit <- NULL; membership <- supplied_membership
    }

    graph_embedding <- if (!is.null(membership_fit) && !is.null(membership_fit$embedding))
      membership_fit$embedding else .dk_graph_spectral_embedding(.dk_as_weighted_graph(graph), dims = wnn_dims)
    graph_cells <- rownames(graph_embedding)
    rna_counts <- Seurat::GetAssayData(object = object, assay = assay, slot = slot)
    atac_counts <- Seurat::GetAssayData(object = object, assay = atac_assay, slot = atac_slot)
    if (!all(graph_cells %in% colnames(rna_counts)) || !all(graph_cells %in% colnames(atac_counts)))
      stop("WNN graph cells are not fully represented in both RNA and ATAC assays", call. = FALSE)
    rna_counts <- rna_counts[, graph_cells, drop = FALSE]
    atac_counts <- atac_counts[, graph_cells, drop = FALSE]
    membership <- .dk_align_membership(membership, graph_cells)
    if (!is.null(membership_fit)) membership_fit$membership <- membership
    group_aligned <- .dk_align_vector(group, graph_cells, "group")
    split_aligned <- .dk_align_vector(split, graph_cells, "split_by")
    hard_stratum <- if (!is.null(group_aligned) || !is.null(split_aligned))
      .dk_stratum(group_aligned, split_aligned, length(graph_cells)) else NULL

    if (is.null(normalize)) normalize <- identical(slot, "counts")
    rna_args <- c(list(x = rna_counts, embedding = graph_embedding,
      membership = if (!is.null(membership_fit)) membership_fit else membership,
      group = group_aligned, split_by = split_aligned, normalize = normalize,
      normalization_scale_factor = normalization_scale_factor), dots)
    rna_res <- do.call(dropout_killer, rna_args)

    atac_res <- .dk_atac_dropout_killer(counts = atac_counts, embedding = graph_embedding,
      membership = membership, membership_fit = membership_fit, group = group_aligned,
      hard_stratum = hard_stratum, bfdr = atac_bfdr, pi_min = atac_pi_min,
      pre_pi = atac_pre_pi, score_min = atac_score_min,
      min_observed_donors = atac_min_observed_donors,
      min_effective_donors = atac_min_effective_donors,
      capture_max_iter = atac_capture_max_iter, capture_tol = atac_capture_tol,
      capture_eps = atac_capture_eps, max_candidates = atac_max_candidates,
      tree_weight = take("tree_weight", 0.5), tree_tau = take("tree_tau", NULL),
      local_k = take("local_k", 30L), candidate_k = take("candidate_k", 100L),
      local_info_kappa = take("local_info_kappa", 5), phi_prior = atac_phi_prior,
      phi_kappa = atac_phi_kappa, phi_floor = atac_phi_floor)

    if (is.null(rna_new_assay)) rna_new_assay <- new_assay
    object[[rna_new_assay]] <- Seurat::CreateAssayObject(data = rna_res$expression)
    object[[atac_new_assay]] <- Seurat::CreateAssayObject(data = atac_res$expression)
    misc <- object@misc
    if (is.null(misc$DropoutKiller)) misc$DropoutKiller <- list()
    misc$DropoutKiller[[rna_new_assay]] <- list(modality = "RNA", settings = rna_res$settings,
      n_memberships = length(unique(membership)), n_detected = nrow(rna_res$events),
      n_recovered = if (nrow(rna_res$events)) sum(rna_res$events$changed) else 0L,
      geometry = "frozen_WNN")
    misc$DropoutKiller[[atac_new_assay]] <- list(modality = "ATAC", settings = atac_res$settings,
      n_memberships = length(unique(membership)), n_candidates = nrow(atac_res$candidate_events),
      n_detected = nrow(atac_res$events),
      n_recovered = if (nrow(atac_res$events)) sum(atac_res$events$changed) else 0L,
      geometry = "frozen_WNN", value_semantics = "expected_accessibility_not_observed_fragments")
    object@misc <- misc

    result <- list(membership = membership, membership_fit = membership_fit,
      embedding = graph_embedding, rna = rna_res, atac = atac_res,
      settings = list(modality = "multiome", rna_assay = assay, atac_assay = atac_assay,
        wnn_graph = wnn_graph, wnn_dims = as.integer(wnn_dims), frozen_geometry = TRUE))
    class(result) <- "DropoutKillerMultiomeResult"
    if (return_result) list(object = object, result = result) else object
  }
}
