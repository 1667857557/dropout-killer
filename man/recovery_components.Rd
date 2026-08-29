\name{masked_factor_prediction}
\alias{masked_factor_prediction}
\alias{weighted_neighbor_prediction}
\alias{recover_dropout_expression}
\alias{sample_dropout_expression}
\title{Selective dropout recovery components}
\description{Predicts supplied dropout-mask events without overwriting observed coordinates. The P1 stabilized-state engine builds target-safe cross-fitted predictors, applies one geometry-aware smoothing step to the predictor state, and performs positive-donor ridge recovery. Tree-local, masked-factor, and Gaussian-neighbor engines remain available for comparison.}
\usage{
masked_factor_prediction(x, membership, mask, factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  cap_quantile = NULL, return_events = FALSE,
  target_mode = c("positive", "all_observed"),
  factor_crossfit_folds = 1L, factor_crossfit_seed = 1L,
  support_adaptive_rank = FALSE, bias_kappa = Inf)

weighted_neighbor_prediction(x, embedding, membership, mask, k = 30L,
  sigma = NULL, min_positive_neighbors = 1L, positive_only = TRUE,
  return_events = FALSE)

recover_dropout_expression(x, mask, membership, embedding = NULL,
  neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, return_details = FALSE,
  recovery_method = c("masked_factor", "p1_stabilized_state",
  "tree_local_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  factor_target = c("positive", "all_observed"),
  membership_fit = NULL, hard_stratum = NULL,
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, factor_crossfit_folds = 1L,
  factor_crossfit_seed = 1L, support_adaptive_rank = FALSE,
  bias_kappa = Inf, predictor_smoothing = 0.25)

sample_dropout_expression(result, n = 1L, seed = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell, or a \code{DropoutKillerMembership} object for direct tree-local recovery.}
\item{mask}{Sparse or dense logical matrix containing only zero-valued dropout events.}
\item{embedding}{Cell-by-dimension biological representation. Required for \code{"p1_stabilized_state"}, \code{"neighbor"}, and \code{"tree_local_factor"}.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor controls. Their positional slots are unchanged.}
\item{cap_quantile}{Optional upper cap used by historical recovery engines.}
\item{return_details}{Return event diagnostics and local geometry with the recovered matrix.}
\item{recovery_method}{Recovery engine. The lower-level API retains \code{"masked_factor"} as its compatibility default. \code{"p1_stabilized_state"} exposes the production high-level engine; \code{"tree_local_factor"} and \code{"neighbor"} retain the comparison estimators.}
\item{factor_rank, factor_features, factor_ridge}{Non-target coexpression-factor and ridge settings.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum target donors required before a factor residual model is attempted. Tree-local recovery otherwise falls back to its query-specific weighted positive mean when local support is adequate.}
\item{target_mode}{Target estimand for \code{masked_factor_prediction}.}
\item{factor_target}{Target estimand for direct/end-to-end recovery. \code{"tree_local_factor"} intentionally supports the positive-conditional target only.}
\item{membership_fit}{Optional \code{DropoutKillerMembership} carrying retained walktrap hierarchies. If omitted, tree-local direct recovery uses embedding-local weighting inside the supplied memberships.}
\item{hard_stratum}{Optional hard biological boundary for direct recovery. Different strata receive exactly zero donor weight.}
\item{tree_weight}{Relative contribution of available walktrap-tree proximity to the local kernel, in [0,1].}
\item{tree_tau}{Positive tree-distance decay scale. \code{NULL} uses a robust within-stratum median positive tree distance.}
\item{local_k}{Neighbor rank defining the adaptive embedding bandwidth for each query cell.}
\item{candidate_k}{Maximum embedding-neighbor component of the candidate pool; the query's final membership is also retained as a high-weight core.}
\item{min_effective_donors}{Minimum Kish effective positive-donor sample size required to recover an event.}
\item{local_info_kappa}{Non-negative shrinkage strength multiplying residual-factor contribution by \eqn{n_eff/(n_eff+kappa)}.}
\item{factor_crossfit_folds}{Number of deterministic target-gene folds. Each target gene is excluded from the factor state used to predict its fold. The legacy default of 1 excludes all current target genes together.}
\item{factor_crossfit_seed}{Integer controlling deterministic target-gene fold assignment.}
\item{support_adaptive_rank}{Logical. If true, cap the target regression rank from its available positive-donor support.}
\item{bias_kappa}{Non-negative LOO bias-shrinkage constant. \code{Inf} disables bias calibration and preserves the legacy result.}
\item{predictor_smoothing}{Weight in \eqn{[0,1]} for the single row-stochastic hierarchy/embedding smoothing step used by \code{"p1_stabilized_state"}.}
\item{return_events}{Return event-level diagnostics from a prediction function.}
\item{k, sigma, positive_only}{Arguments for \code{weighted_neighbor_prediction}.}
\item{result}{A \code{DropoutKillerResult} with predictive uncertainty.}
\item{n}{Number of uncertainty-aware completed matrices.}
\item{seed}{Optional random seed.}
}
\details{For P1 stabilized-state recovery, each target-gene fold is excluded before predictor construction. The standardized non-target predictor matrix receives exactly one row-stochastic smoothing step over the hierarchy/embedding geometry before SVD, while the final expression values are never smoothed. Positive donors then train the P1 ridge state with analytic leave-one-out shrinkage and optional bias calibration. For tree-local recovery, hard strata define absolute borrowing boundaries and nearby sibling memberships may contribute through a decaying tree/embedding kernel. Positive completed-matrix draws use Gamma moment matching.}
\value{\code{recover_dropout_expression} returns a selectively recovered matrix or, with \code{return_details=TRUE}, a list containing expression, events, uncertainty status, and local geometry for the tree-local engine. Event diagnostics include local positive mean/variance/prevalence, effective donors, tree and embedding distances, factor predictability, shrinkage, and prediction standard deviation.}
