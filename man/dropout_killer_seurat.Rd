\name{dropout_killer_seurat}
\alias{dropout_killer_seurat}
\title{Run DropoutKiller on a Seurat object}
\description{Runs Lean hierarchy detection and selective recovery on a Seurat object. Paired RNA/ATAC inputs automatically recompute PCA, LSI excluding component one, and WNN.}
\usage{
dropout_killer_seurat(object, assay = NULL, slot = "counts",
  reduction = "pca", dims = 1:20, group_by = NULL, split_by = NULL,
  new_assay = "DropoutKiller", return_result = FALSE,
  normalize = NULL, normalization_scale_factor = 1e4,
  modality = c("auto", "rna", "wnn"), atac_assay = NULL,
  lean_model = NULL, lean_control = list(), wnn_npcs = 40L, wnn_k = 30L, ...)
}
\arguments{
\item{modality}{Auto detects paired RNA/ATAC, rna opts out, wnn requires ATAC.}
\item{atac_assay}{Explicit ATAC assay name if detection is ambiguous.}
\item{lean_model, lean_control}{Fitted model and calibration options.}
\item{wnn_npcs}{Maximum recomputed PCA/LSI dimensions. LSI1 is always excluded.}
\item{wnn_k}{WNN neighbors within broad classes.}
\item{object}{Seurat object.}
\item{assay}{Input assay; defaults to RNA when available, otherwise the active assay.}
\item{slot}{Input slot. Defaults to \code{"counts"}; when \code{"data"} or another pre-normalized slot is selected, automatic normalization defaults to false.}
\item{reduction}{Low-dimensional reduction name.}
\item{dims}{Dimensions to use.}
\item{group_by}{Metadata column for broad biological classes, not fine subclusters. Required for automatic proxy-negative calibration.}
\item{split_by}{Optional metadata column for an additional hard recovery boundary. Rejected by the new Lean route, which uses broad classes only.}
\item{new_assay}{Name of the recovered assay.}
\item{return_result}{Return both the modified object and the full DropoutKiller result.}
\item{normalize}{Optional logical. \code{NULL} chooses \code{TRUE} for \code{slot="counts"} and \code{FALSE} otherwise.}
\item{normalization_scale_factor}{Positive target library size before \code{log1p}; default \code{1e4}, matching ALRA.}
\item{...}{Additional arguments passed to \code{dropout_killer}.}
}
\value{A Seurat object, or a list containing the object and result when \code{return_result=TRUE}. Recovered values are stored on the normalized continuous expression scale, not as integer counts.}
