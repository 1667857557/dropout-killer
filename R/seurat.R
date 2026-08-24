#' Run DropoutKiller on a Seurat object
#'
#' The recovered matrix is written to the data slot of a new assay. It is not
#' presented as raw counts because weighted recovery produces continuous values.
#'
#' @export
dropout_killer_seurat <- function(object, assay = NULL, slot = "data", reduction = "pca",
                                  dims = 1:20, group_by = NULL, split_by = NULL,
                                  new_assay = "DropoutKiller", return_result = FALSE, ...) {
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
  if (identical(slot, "counts")) warning("input slot='counts': recovered values are continuous and will still be stored as assay data, not counts", call. = FALSE)
  res <- dropout_killer(x = x, embedding = emb, group = group, split_by = split, ...)
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
