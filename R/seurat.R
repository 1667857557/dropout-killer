#' Run DropoutKiller on a Seurat object
#'
#' By default the raw `counts` slot is supplied to `dropout_killer()`, which
#' applies the ALRA library-size normalization to 10,000 counts per cell followed
#' by `log1p`. The recovered matrix is written to the data slot of a new assay;
#' it is not presented as raw counts because recovery produces continuous values.
#'
#' If a pre-normalized slot such as `data` is selected explicitly, normalization
#' defaults to `FALSE` to avoid applying the ALRA transform twice. This can be
#' overridden with `normalize`.
#'
#' @export
dropout_killer_seurat <- function(object, assay = NULL, slot = "counts", reduction = "pca",
                                  dims = 1:20, group_by = NULL, split_by = NULL,
                                  new_assay = "DropoutKiller", return_result = FALSE,
                                  normalize = NULL, normalization_scale_factor = 1e4, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required", call. = FALSE)
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object", call. = FALSE)
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  x <- Seurat::GetAssayData(object = object, assay = assay, slot = slot)
  emb <- Seurat::Embeddings(object, reduction = reduction)
  dims <- as.integer(dims); dims <- dims[dims >= 1L & dims <= ncol(emb)]
  if (!length(dims)) stop("no requested dimensions exist in reduction", call. = FALSE)
  emb <- emb[, dims, drop = FALSE]
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
  if (is.null(normalize)) normalize <- identical(slot, "counts")
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) stop("normalize must be TRUE or FALSE", call. = FALSE)
  if (identical(slot, "counts") && !normalize)
    warning("slot='counts' with normalize=FALSE leaves raw counts on the recovery scale", call. = FALSE)
  if (!identical(slot, "counts") && normalize)
    warning("normalizing a non-count Seurat slot; set normalize=FALSE if the selected slot is already library/log normalized", call. = FALSE)
  res <- dropout_killer(
    x = x, embedding = emb, group = group, split_by = split,
    normalize = normalize, normalization_scale_factor = normalization_scale_factor, ...
  )
  object[[new_assay]] <- Seurat::CreateAssayObject(data = res$expression)
  misc <- object@misc
  if (is.null(misc$DropoutKiller)) misc$DropoutKiller <- list()
  misc$DropoutKiller[[new_assay]] <- list(settings = res$settings,
                                           n_memberships = length(unique(res$membership)),
                                           n_detected = nrow(res$events),
                                           n_recovered = if (nrow(res$events)) sum(res$events$changed) else 0L)
  object@misc <- misc
  if (return_result) list(object = object, result = res) else object
}
