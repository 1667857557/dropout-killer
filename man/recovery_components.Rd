\name{masked_factor_prediction}
\alias{masked_factor_prediction}
\alias{weighted_neighbor_prediction}
\alias{recover_dropout_expression}
\alias{sample_dropout_expression}
\title{Selective dropout recovery components}
\description{Predicts supplied dropout-mask events without overwriting observed coordinates. The tree-local engine combines retained SuperCell hierarchy proximity with the original biological embedding to build query-specific positive-expression baselines, then predicts only residual expression with locally weighted non-target coexpression factors. The historical masked-factor and Gaussian-neighbor engines remain available for comparison.}
\usage{
masked_factor_prediction(x, membership, mask, factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  cap_quantile = NULL, return_events = FALSE,
  target_mode = c("positive", "all_observed"))

weighted_neighbor_prediction(x, embedding, membership, mask, k = 30L,
  sigma = NULL, min_positive_neighbors = 1L, positive_only = TRUE,
  return_events = FALSE)

recover_dropout_expression(x, mask, membership, embedding = NULL,
  neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, return_details = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  factor_target = c("positive", "all_observed"),
  membership_fit = NULL, hard_stratum = NULL,
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5)

sample_dropout_expression(result, n = 1L, seed = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell, or a \code{DropoutKillerMembership} object for direct tree-local recovery.}
\item{mask}{Sparse or dense logical matrix containing only zero-valued dropout events.}
\item{embedding}{Cell-by-dimension biological representation. Required for \code{"neighbor"} and \code{"tree_local_factor"}.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor controls. Their positional slots are unchanged.}
\item{cap_quantile}{Optional upper cap used by historical recovery engines.}
\item{return_details}{Return event diagnostics and local geometry with the recovered matrix.}
\item{recovery_method}{Recovery engine. \code{"masked_factor"} remains the default comparator; \code{"tree_local_factor"} uses hierarchy- and embedding-weighted local recovery; \code{"neighbor"} uses the legacy Gaussian cell-neighbor estimator.}
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
\item{return_events}{Return event-level diagnostics from a prediction function.}
\item{k, sigma, positive_only}{Arguments for \code{weighted_neighbor_prediction}.}
\item{result}{A \code{DropoutKillerResult} with predictive uncertainty.}
\item{n}{Number of uncertainty-aware completed matrices.}
\item{seed}{Optional random seed.}
}
\details{For tree-local recovery, hard strata define absolute borrowing boundaries. A walktrap tree cut is not treated as an absolute wall: nearby sibling memberships may contribute through a decaying tree/embedding kernel. The local positive baseline for target gene g and query cell c is \eqn{sum_j w_cj X_gj / sum_j w_cj} over reliable positive donors. Kish effective donor size is \eqn{(sum w)^2/sum(w^2)}. Non-target factors predict leave-one-out residuals beyond this local baseline using query-weighted ridge regression; unsupported residual structure collapses to the query-specific local mean, not the whole-membership mean. Predictive uncertainty includes held-out residual error, local-mean uncertainty, and local factor leverage. Positive completed-matrix draws use Gamma moment matching.}
\value{\code{recover_dropout_expression} returns a selectively recovered matrix or, with \code{return_details=TRUE}, a list containing expression, events, uncertainty status, and local geometry for the tree-local engine. Event diagnostics include local positive mean/variance/prevalence, effective donors, tree and embedding distances, factor predictability, shrinkage, and prediction standard deviation.}
