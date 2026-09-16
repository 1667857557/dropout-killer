\name{dropout_killer_seurat}
\alias{dropout_killer_seurat}
\title{Run DropoutKiller on a Seurat object}
\description{Runs the current Lean hierarchy detector and selective recovery on a Seurat object. Paired RNA+ATAC inputs automatically rebuild RNA PCA, ATAC LSI excluding component one, and WNN within broad classes.}
\usage{
dropout_killer_seurat(object, assay = NULL, slot = "counts",
  group_by = NULL, new_assay = "DropoutKiller",
  return_result = FALSE, normalize = NULL,
  normalization_scale_factor = 1e4,
  modality = c("auto", "rna", "wnn"), atac_assay = NULL,
  lean_model = NULL, lean_control = list(), wnn_npcs = 40L,
  wnn_k = 30L, ...)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{RNA input assay; defaults to RNA when available.}
\item{slot}{Input slot. Defaults to \code{"counts"}.}
\item{group_by}{Metadata column containing broad biological classes, not fine subclusters. Required for automatic proxy-negative calibration.}
\item{new_assay}{Name of the recovered assay.}
\item{return_result}{Return both the modified object and the full DropoutKiller result.}
\item{normalize}{Optional logical. NULL chooses TRUE for the counts slot and FALSE otherwise.}
\item{normalization_scale_factor}{Target library size before \code{log1p}; the Lean contract requires 10000.}
\item{modality}{\code{"auto"} uses paired RNA+ATAC when unambiguous, \code{"rna"} explicitly opts out of ATAC, and \code{"wnn"} requires ATAC.}
\item{atac_assay}{Explicit ATAC assay name when more than one candidate exists.}
\item{lean_model}{Optional pre-fitted Lean detector model of the matching RNA/WNN modality.}
\item{lean_control}{Named calibration arguments.}
\item{wnn_npcs}{Maximum recomputed RNA PCA / ATAC LSI dimensions. LSI1 is always excluded.}
\item{wnn_k}{WNN neighbor count within each broad class.}
\item{...}{Additional current-workflow arguments passed to \code{dropout_killer}.}
}
\value{A Seurat object, or a list containing the object and full result when \code{return_result=TRUE}. Recovered values are stored on the normalized continuous expression scale in the new assay's data slot.}
