\name{dropout_killer_seurat}
\alias{dropout_killer_seurat}
\title{Run DropoutKiller on a Seurat object}
\description{Uses a Seurat assay and reduction, runs the selective recovery workflow, and writes continuous recovered expression to a new assay data slot.}
\usage{
dropout_killer_seurat(object, assay = NULL, slot = "data",
  reduction = "pca", dims = 1:20, group_by = NULL, split_by = NULL,
  new_assay = "DropoutKiller", return_result = FALSE, ...)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{Input assay; defaults to the active assay.}
\item{slot}{Input slot, normally \code{"data"}.}
\item{reduction}{Low-dimensional reduction name.}
\item{dims}{Dimensions to use.}
\item{group_by}{Optional metadata column for hard major-cell-type strata.}
\item{split_by}{Optional metadata column for an additional hard boundary.}
\item{new_assay}{Name of the recovered assay.}
\item{return_result}{Return both the modified object and the full DropoutKiller result.}
\item{...}{Additional arguments passed to \code{dropout_killer}.}
}
\value{A Seurat object, or a list containing the object and result when \code{return_result=TRUE}.}
