#' Run DropoutKiller on a Seurat object
#'
#' Raw RNA counts are used by default and recovered values are written to the
#' data slot of a new assay. RNA-only inputs rebuild detector geometry from the
#' normalized RNA SVD. Paired RNA+ATAC inputs rebuild RNA PCA, ATAC TF-IDF/LSI
#' (excluding LSI1), and WNN within each supplied broad class.
#'
#' @export
dropout_killer_seurat <- function(object, assay = NULL, slot = "counts",
                                  group_by = NULL, new_assay = "DropoutKiller",
                                  return_result = FALSE, normalize = NULL,
                                  normalization_scale_factor = 1e4,
                                  modality = c("auto", "rna", "wnn"),
                                  atac_assay = NULL, lean_model = NULL,
                                  lean_control = list(), wnn_npcs = 40L,
                                  wnn_k = 30L, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required", call. = FALSE)
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object", call. = FALSE)
  modality <- match.arg(modality)
  dots <- list(...)
  selected_atac <- if (modality != "rna") .dk_atac_assay(object, atac_assay) else NULL
  use_wnn <- modality == "wnn" || (modality == "auto" && !is.null(selected_atac))
  if (is.null(assay)) assay <- if ("RNA" %in% names(object@assays)) "RNA" else SeuratObject::DefaultAssay(object)
  if (use_wnn && identical(assay, selected_atac)) stop("select the RNA assay for dropout detection", call. = FALSE)
  expected_method <- .dk_lean_method(use_wnn)
  if (!is.null(dots$detection_method) && !identical(dots$detection_method, expected_method))
    stop("explicit detector conflicts with input modality; set modality explicitly", call. = FALSE)
  dots$detection_method <- expected_method

  x <- .dk_seurat_layer(object, assay, slot)
  meta <- object[[]]
  group <- NULL
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(meta)) stop("group_by column not found", call. = FALSE)
    group <- meta[[group_by]]; names(group) <- rownames(meta)
  }
  if (is.null(normalize)) normalize <- identical(slot, "counts")
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) stop("normalize must be TRUE or FALSE", call. = FALSE)
  if (identical(slot, "counts") && !normalize)
    warning("slot='counts' with normalize=FALSE leaves raw counts on the recovery scale", call. = FALSE)
  if (!identical(slot, "counts") && normalize)
    warning("normalizing a non-count Seurat slot; set normalize=FALSE if the selected slot is already library/log normalized", call. = FALSE)

  value <- function(name, default) if (is.null(dots[[name]])) default else dots[[name]]
  geometry <- NULL
  if (use_wnn) {
    build <- function(obj) build_wnn_supercell(
      obj, rna_assay = assay, atac_assay = selected_atac,
      group = group, gamma = value("gamma", 150), npcs = wnn_npcs,
      k_nn = wnn_k, seed = value("seed", 12345L)
    )
    geometry <- build(object)
    emb <- geometry$embedding
    if (is.null(lean_model)) {
      if (!identical(slot, "counts") || !normalize) stop("automatic WNN calibration requires raw RNA counts", call. = FALSE)
      rebuild <- function(counts) {
        obj <- object
        obj[[assay]] <- SeuratObject::CreateAssayObject(counts = counts)
        build(obj)
      }
      lean_model <- do.call(calibrate_lean_detector, c(list(
        counts = x, group = group, gamma = value("gamma", 150),
        neighbor_k = value("neighbor_k", 30L), rank = value("rank", "auto"),
        geometry_builder = rebuild
      ), lean_control))
    }
  } else {
    y <- if (normalize) .dk_alra_library_log(x) else x
    emb <- .dk_lean_svd(y, value("rank", "auto"), value("seed", 12345L))$embedding
  }

  res <- do.call(dropout_killer, c(list(
    x = x, embedding = emb, group = group,
    normalize = normalize, normalization_scale_factor = normalization_scale_factor,
    lean_model = lean_model, lean_control = lean_control, lean_geometry = geometry
  ), dots))
  if (!is.null(geometry)) res$settings$wnn <- geometry$provenance
  object[[new_assay]] <- Seurat::CreateAssayObject(data = res$expression)
  misc <- object@misc
  if (is.null(misc$DropoutKiller)) misc$DropoutKiller <- list()
  misc$DropoutKiller[[new_assay]] <- list(
    settings = res$settings,
    n_memberships = length(unique(res$membership)),
    n_detected = nrow(res$events),
    n_recovered = if (nrow(res$events)) sum(res$events$changed) else 0L
  )
  object@misc <- misc
  if (return_result) list(object = object, result = res) else object
}
