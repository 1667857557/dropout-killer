\name{masked_factor_prediction}
\alias{masked_factor_prediction}
\alias{weighted_neighbor_prediction}
\alias{recover_dropout_expression}
\alias{sample_dropout_expression}
\title{Selective dropout recovery components}
\description{Predicts supplied dropout-mask events without overwriting observed coordinates. The barycentric engine treats final membership as a soft biological prior inside an absolute hard stratum and learns target-independent simplex donor weights from the biological embedding, while the tree-local engine uses membership-local hierarchy/embedding weighting plus residual factors. Historical masked-factor and Gaussian-neighbor engines remain available for comparison.}
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
  recovery_method = c("masked_factor", "tree_local_factor", "barycentric", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  factor_target = c("positive", "all_observed"),
  membership_fit = NULL, hard_stratum = NULL,
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, membership_penalty = 1,
  barycentric_lambda = 1, barycentric_iter = 20L,
  barycentric_tol = 1e-6)

sample_dropout_expression(result, n = 1L, seed = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell, or a \code{DropoutKillerMembership} object carrying the retained hierarchy.}
\item{mask}{Sparse or dense logical matrix containing only zero-valued dropout events.}
\item{embedding}{Cell-by-dimension biological representation. Required for \code{"neighbor"}, \code{"tree_local_factor"}, and \code{"barycentric"}.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor controls. Their positional slots are unchanged.}
\item{cap_quantile}{Optional upper cap used by historical recovery engines.}
\item{return_details}{Return event diagnostics and local geometry with the recovered matrix.}
\item{recovery_method}{Recovery engine. \code{"masked_factor"} remains the default comparator; \code{"tree_local_factor"} uses hierarchy- and embedding-weighted recovery inside final memberships; \code{"barycentric"} learns query-cell simplex donor weights across local embedding candidates inside hard strata with membership/hierarchy soft priors; \code{"neighbor"} uses the legacy Gaussian estimator.}
\item{factor_rank, factor_features, factor_ridge}{Non-target coexpression-factor and ridge settings for factor engines.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum target donors required before a factor residual model is attempted.}
\item{target_mode}{Target estimand for \code{masked_factor_prediction}.}
\item{factor_target}{Target estimand for direct/end-to-end recovery. \code{"tree_local_factor"} and \code{"barycentric"} intentionally support the positive-conditional target only.}
\item{membership_fit}{Optional \code{DropoutKillerMembership} carrying retained walktrap hierarchies. If omitted, barycentric recovery still uses embedding and membership priors but lacks hierarchy distance.}
\item{hard_stratum}{Optional hard biological boundary for direct recovery. Different strata receive exactly zero donor weight under tree-local and barycentric recovery.}
\item{tree_weight}{Relative contribution of available walktrap-tree proximity to the local geometry prior, in [0,1].}
\item{tree_tau}{Positive tree-distance decay scale. \code{NULL} uses a robust median positive distance where hierarchy information is available.}
\item{local_k}{Neighbor rank defining the adaptive embedding bandwidth for each query cell.}
\item{candidate_k}{Maximum embedding-neighbor candidate-pool size. Barycentric recovery may include candidates from nearby sibling memberships within the same hard stratum.}
\item{min_effective_donors}{Minimum Kish effective positive-donor sample size required to recover an event.}
\item{local_info_kappa}{Non-negative shrinkage strength used by the tree-local residual factor.}
\item{membership_penalty}{Non-negative barycentric donor cost for crossing the final membership cut. It is a preference rather than a hard boundary.}
\item{barycentric_lambda}{Positive strength shrinking learned barycentric weights toward the hierarchy/embedding/membership prior.}
\item{barycentric_iter}{Maximum projected-gradient iterations for each query-cell simplex problem.}
\item{barycentric_tol}{Positive convergence tolerance for barycentric weight updates.}
\item{return_events}{Return event-level diagnostics from a prediction function.}
\item{k, sigma, positive_only}{Arguments for \code{weighted_neighbor_prediction}.}
\item{result}{A \code{DropoutKillerResult} with predictive uncertainty.}
\item{n}{Number of uncertainty-aware completed matrices.}
\item{seed}{Optional random seed.}
}
\details{For barycentric recovery, explicit hard strata define absolute borrowing boundaries. Inside a hard stratum, embedding k-nearest neighbors define a small donor candidate set. A geometry prior combines embedding distance, available hierarchy distance, and a finite penalty for crossing the final membership cut. Learned donor weights solve a non-negative, unit-sum convex reconstruction problem in embedding space and are independent of target-gene values. For a dropout target g in query cell c, these state weights are renormalized over donors with reliable positive X_gj and used to estimate the positive-conditional expression. Kish effective donor size is \eqn{1 / sum_j w_j^2} after target-specific renormalization. Weighted local variance supplies a working predictive standard deviation. Unsupported events remain unchanged rather than being forced to a value.}
\value{\code{recover_dropout_expression} returns a selectively recovered matrix or, with \code{return_details=TRUE}, a list containing expression, events, uncertainty status, and engine-specific local geometry. Barycentric geometry stores a sparse query-cell-by-donor weight matrix together with the prior matrix, bandwidth, reconstruction error, cross-membership weight summary, and hierarchy/embedding distance diagnostics.}
