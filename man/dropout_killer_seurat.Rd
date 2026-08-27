\name{dropout_killer_seurat}
\alias{dropout_killer_seurat}
\title{Run DropoutKiller on a Seurat object}
\description{Uses a Seurat assay and reduction, runs the selective recovery workflow, and writes continuous recovered expression to a new assay data slot. The raw counts slot is used by default and is internally transformed with ALRA-style library-size normalization and \code{log1p}. Under the default detector, \code{group_by} defines the major cell-class blocks for global ALRA zero detection; \code{split_by} remains a recovery boundary.}
\usage{
dropout_killer_seurat(object, assay = NULL, slot = "counts",
  reduction = "pca", dims = 1:20, group_by = NULL, split_by = NULL,
  new_assay = "DropoutKiller", return_result = FALSE,
  normalize = NULL, normalization_scale_factor = 1e4, ...)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{Input assay; defaults to the active assay.}
\item{slot}{Input slot. Defaults to \code{"counts"}; when \code{"data"} or another pre-normalized slot is selected, automatic normalization defaults to false.}
\item{reduction}{Low-dimensional reduction name.}
\item{dims}{Dimensions to use.}
\item{group_by}{Optional metadata column identifying the major cell class. With the default \code{"alra_global_by_group"} detector, each value is one independent global ALRA zero-detection block. If omitted, all cells form one detection block. It also remains a hard recovery stratum.}
\item{split_by}{Optional metadata column for an additional hard recovery boundary. It does not split default ALRA detection blocks.}
\item{new_assay}{Name of the recovered assay.}
\item{return_result}{Return both the modified object and the full DropoutKiller result.}
\item{normalize}{Optional logical. \code{NULL} chooses \code{TRUE} for \code{slot="counts"} and \code{FALSE} otherwise.}
\item{normalization_scale_factor}{Positive target library size before \code{log1p}; default \code{1e4}, matching ALRA.}
\item{...}{Additional arguments passed to \code{dropout_killer}.}
}
\value{A Seurat object, or a list containing the object and result when \code{return_result=TRUE}. Recovered values are stored on the normalized continuous expression scale, not as integer counts.}
