\name{dropout_killer_seurat}
\alias{dropout_killer_seurat}
\title{Run DropoutKiller on a Seurat object}
\description{
Uses a Seurat assay and reduction, runs the selective recovery workflow, and
writes continuous recovered expression to a new assay data slot. The default
scGACL detector requires raw counts and a \code{group_by} cell-subpopulation
column. Detection applies \code{log(1.01 + raw_count)} internally, while the
recovery working matrix is independently ALRA library-size normalized and
\code{log1p} transformed when \code{normalize=TRUE}.
}
\usage{
dropout_killer_seurat(object, assay = NULL, slot = "counts",
  reduction = "pca", dims = 1:20, group_by = NULL, split_by = NULL,
  new_assay = "DropoutKiller", return_result = FALSE,
  normalize = NULL, normalization_scale_factor = 1e4, ...)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{Input assay; defaults to the active assay.}
\item{slot}{Input slot. Defaults to \code{"counts"}. A pre-normalized \code{"data"} slot is not a source-faithful input for the default scGACL detector; explicitly choose a compatible detector when using such a slot.}
\item{reduction}{Low-dimensional reduction name.}
\item{dims}{Dimensions to use.}
\item{group_by}{Metadata column defining cell subpopulations for the default scGACL detector, typically major cell type. It also remains a hard recovery stratum.}
\item{split_by}{Optional metadata column for an additional hard recovery boundary.}
\item{new_assay}{Name of the recovered assay.}
\item{return_result}{Return both the modified object and the full DropoutKiller result.}
\item{normalize}{Optional logical. \code{NULL} chooses \code{TRUE} for \code{slot="counts"} and \code{FALSE} otherwise. This controls the recovery working scale, not the raw-count scGACL detector transform.}
\item{normalization_scale_factor}{Positive target library size before recovery \code{log1p}; default \code{1e4}.}
\item{...}{Additional arguments passed to \code{dropout_killer}.}
}
\value{A Seurat object, or a list containing the object and result when \code{return_result=TRUE}. Recovered values are stored on the normalized continuous expression scale, not as integer counts.}
