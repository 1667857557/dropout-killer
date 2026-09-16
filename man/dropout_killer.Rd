\name{dropout_killer}
\alias{dropout_killer}
\title{Run the current DropoutKiller workflow}
\description{Detects recoverable zeros with the calibrated SuperCell hierarchy Lean detector and selectively reconstructs selected events. RNA+ATAC WNN geometry is supplied by the Seurat wrapper.}
\usage{
dropout_killer(x, embedding, group = NULL, gamma = 150,
  k_knn = 5L, rank = "auto", neighbor_k = 30L,
  neighbor_sigma = NULL, min_positive_neighbors = 1L,
  neighbor_positive_only = TRUE, cap_quantile = NULL,
  seed = 12345L, return_score = FALSE,
  recovery_method = c("p1_stabilized_state", "masked_factor",
  "tree_local_factor", "neighbor"), factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 2,
  min_feature_observed = 20L, min_target_observed = 8L,
  detection_method = c("Supercell_hierarchy_Lean_membership",
  "Supercell_hierarchy_Lean_membership_WNN"),
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, factor_crossfit_folds = 5L,
  factor_crossfit_seed = 1L, support_adaptive_rank = TRUE,
  bias_kappa = 10, predictor_smoothing = 0.25,
  lean_model = NULL, lean_control = list(), lean_geometry = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell RNA count/expression matrix. Raw counts are expected when automatic calibration is used.}
\item{embedding}{Cell-by-dimension representation used by the recovery engine. RNA-only detection rebuilds its own normalized-RNA SVD geometry.}
\item{group}{Broad biological class per cell. It is the hard boundary for hierarchy construction and borrowing. Automatic calibration requires at least two suitable broad classes.}
\item{gamma}{SuperCell graining level; the Lean hierarchy uses a floor \code{n/gamma} cut subject to connected components.}
\item{k_knn}{RNA kNN graph neighbor count.}
\item{rank}{Low-rank SVD rank or \code{"auto"}.}
\item{neighbor_k}{Maximum hierarchy support-neighbor count.}
\item{neighbor_sigma, min_positive_neighbors, neighbor_positive_only, cap_quantile}{Controls used by the optional comparison recovery engines.}
\item{seed}{Random seed for deterministic stochastic components.}
\item{return_score}{Also return a sparse score matrix over called zero events.}
\item{recovery_method}{Recovery engine. \code{"p1_stabilized_state"} is the production default; the remaining choices are retained for controlled recovery ablation.}
\item{factor_rank, factor_features, factor_ridge}{Target-safe coexpression-state and ridge settings.}
\item{min_feature_observed}{Minimum observations for a predictor feature.}
\item{min_target_observed}{Minimum positive target donors before a factor state is attempted.}
\item{detection_method}{Lean detector modality. RNA is the default; WNN requires matching precomputed WNN geometry, normally built by \code{dropout_killer_seurat}.}
\item{factor_target}{Recovery target. P1 stabilized state and tree-local recovery support the positive-conditional target.}
\item{tree_weight, tree_tau, local_k, candidate_k, min_effective_donors, local_info_kappa}{Expert controls for hierarchy-aware comparison recovery.}
\item{normalize}{Apply library-size normalization followed by \code{log1p}. Automatic calibration requires \code{TRUE}.}
\item{normalization_scale_factor}{Target library size. The Lean feature contract requires 10000.}
\item{factor_crossfit_folds}{Number of deterministic target-gene folds used to prevent target leakage.}
\item{factor_crossfit_seed}{Seed controlling target-gene fold assignment.}
\item{support_adaptive_rank}{Cap target regression rank using available positive-donor support.}
\item{bias_kappa}{Leave-one-out bias-shrinkage constant.}
\item{predictor_smoothing}{Weight in \eqn{[0,1]} for one row-stochastic geometry smoothing step applied to the target-safe predictor state.}
\item{lean_model}{A fitted \code{DropoutKillerLeanModel}. If NULL, the RNA route calibrates from raw counts.}
\item{lean_control}{Named arguments passed to \code{calibrate_lean_detector}.}
\item{lean_geometry}{Optional precomputed Lean geometry. WNN geometry must match the detector model modality.}
}
\details{Default preprocessing is \deqn{Y_{gc}=\log(1+10^4 X_{gc}/\sum_hX_{hc}).} The Lean detector combines standardized ALRA low-rank evidence, membership-first hierarchy neighbor support, and shrunk membership-versus-group context with a class-balanced logistic model. A zero is called only when its linear score exceeds the calibration threshold; its sigmoid is not interpreted as a biological dropout posterior. P1 recovery excludes target-gene folds from predictor construction, applies one geometry-aware smoothing step, and fits positive-donor ridge states with analytic leave-one-out shrinkage and bias calibration. Observed nonzero coordinates remain unchanged.}
\value{A \code{DropoutKillerResult} containing recovered normalized expression, sparse mask, event table, memberships and retained hierarchy, detection model/provenance, recovery settings, and predictive variance when available.}
