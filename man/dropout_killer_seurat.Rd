\name{dropout_killer_seurat}
\alias{dropout_killer_seurat}
\title{Run DropoutKiller on a Seurat or 10x Multiome object}
\description{Runs the selective recovery workflow on a Seurat object. The default \code{modality="rna"} path preserves the historical RNA behavior. With \code{modality="multiome"}, one frozen WNN-derived membership is shared by RNA and ATAC: RNA uses the existing group-wise ALRA detector and P1 recovery, while ATAC uses a latent accessibility/capture model and selective local Gamma-Poisson empirical-Bayes recovery. Raw assays and fragment evidence are never overwritten.}
\usage{
dropout_killer_seurat(object, assay = NULL, slot = "counts",
  reduction = "pca", dims = 1:20, group_by = NULL, split_by = NULL,
  new_assay = "DropoutKiller", return_result = FALSE,
  normalize = NULL, normalization_scale_factor = 1e4,
  modality = c("rna", "multiome"), atac_assay = "ATAC",
  atac_slot = "counts", wnn_graph = "wsnn", wnn_dims = 30L,
  rna_new_assay = NULL, atac_new_assay = "ATAC_DK",
  atac_bfdr = 0.01, atac_pi_min = 0.8, atac_pre_pi = 0.5,
  atac_score_min = 0.5, atac_min_observed_donors = 2L,
  atac_min_effective_donors = 5, atac_capture_max_iter = 25L,
  atac_capture_tol = 1e-5, atac_capture_eps = 1e-4,
  atac_max_candidates = 5e6, atac_phi_prior = 1,
  atac_phi_kappa = 10, atac_phi_floor = 1e-4, ...)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{RNA/input assay; defaults to the active assay. In multiome mode this is the RNA assay.}
\item{slot}{RNA input slot.}
\item{reduction}{Low-dimensional reduction name used only in RNA-only mode.}
\item{dims}{Dimensions used only in RNA-only mode.}
\item{group_by}{Optional metadata column defining major biological classes. It defines RNA ALRA blocks, ATAC capture-estimation blocks, and a hard recovery stratum.}
\item{split_by}{Optional metadata column for an additional hard recovery boundary.}
\item{new_assay}{Historical RNA recovered-assay name. In multiome mode it is used for RNA unless \code{rna_new_assay} is supplied.}
\item{return_result}{Return both the modified object and full result.}
\item{normalize}{RNA normalization control.}
\item{normalization_scale_factor}{RNA target library size before \code{log1p}.}
\item{modality}{Either \code{"rna"} or \code{"multiome"}.}
\item{atac_assay}{ATAC assay name in multiome mode.}
\item{atac_slot}{ATAC quantitative count slot.}
\item{wnn_graph}{Seurat WNN graph name, typically \code{"wsnn"}.}
\item{wnn_dims}{Spectral dimensions retained from the WNN graph for local recovery geometry.}
\item{rna_new_assay}{Optional RNA recovered-assay name in multiome mode.}
\item{atac_new_assay}{ATAC recovered-assay name. Values are expected accessibility, not observed fragments.}
\item{atac_bfdr}{Bayesian FDR target for high-confidence ATAC dropout selection.}
\item{atac_pi_min}{Minimum leave-one-out membership accessibility probability for an ATAC zero candidate.}
\item{atac_pre_pi}{Membership accessibility prefilter before event construction.}
\item{atac_score_min}{Minimum posterior dropout probability before BFDR ranking.}
\item{atac_min_observed_donors}{Minimum observed nonzero donor cells for a candidate peak in a membership.}
\item{atac_min_effective_donors}{Minimum effective donor support for ATAC detection/recovery.}
\item{atac_capture_max_iter}{Maximum PIC-style coordinate-descent iterations for cell capture probabilities.}
\item{atac_capture_tol}{Convergence tolerance for ATAC capture estimation.}
\item{atac_capture_eps}{Numerical probability boundary.}
\item{atac_max_candidates}{Safety limit on retained ATAC candidate zero events.}
\item{atac_phi_prior}{Prior squared coefficient of variation for local Gamma-Poisson accessibility rates.}
\item{atac_phi_kappa}{Support shrinkage for peak-specific \code{phi} estimates.}
\item{atac_phi_floor}{Minimum peak-specific \code{phi}.}
\item{...}{Additional RNA DropoutKiller arguments. Membership/hierarchy controls such as \code{gamma}, \code{tree_weight}, \code{local_k}, and \code{candidate_k} are also reused by the ATAC branch.}
}
\value{In RNA-only mode, a Seurat object or object/result list as before. In multiome mode, the object receives separate RNA and ATAC recovered assays; when \code{return_result=TRUE}, the result contains the shared frozen membership plus separate \code{rna} and \code{atac} components. ATAC recovered values are continuous expected accessibility values and are not raw fragment counts.}
