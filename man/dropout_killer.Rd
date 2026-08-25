\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Builds or accepts biologically restricted cell memberships, detects selected technical-zero candidates, and selectively recovers only those events. Version 0.6 can retain the full SuperCell walktrap hierarchy and combine tree proximity with the original biological embedding so closer positive donors contribute more strongly to query-specific recovery. Historical masked-factor and neighbor engines remain available. Observed values are never overwritten.}
\usage{
dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5)

DropoutKiller(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5)

run_dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix, typically normalized/log-transformed expression rather than raw counts.}
\item{embedding}{Cell-by-dimension biological representation used for membership construction and local recovery geometry.}
\item{membership}{Optional membership vector or \code{DropoutKillerMembership}. Internally built memberships retain the full walktrap hierarchy.}
\item{group}{Optional hard biological stratum such as major cell type. Different groups never borrow under tree-local recovery.}
\item{split_by}{Optional additional hard boundary such as condition or donor.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters; \code{gamma} controls the final tree cut but the retained hierarchy remains available to recovery.}
\item{rank, max_rank, rank_z, quantile_prob, threshold, min_cells, min_negative}{Dropout-detection controls.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor controls; positional slots are retained.}
\item{cap_quantile}{Optional historical target-gene cap.}
\item{seed}{Random seed.}
\item{return_score}{Also return the sparse detection-score matrix.}
\item{recovery_method}{Recovery engine. \code{"masked_factor"} remains the default comparator. \code{"tree_local_factor"} uses retained hierarchy plus embedding-local positive baselines and locally weighted residual factors. \code{"neighbor"} is the legacy cell-neighbor engine.}
\item{factor_rank, factor_features, factor_ridge}{Non-target coexpression-factor and ridge parameters.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum positive target donors before factor residual learning is attempted.}
\item{detection_method}{Detection engine. Default \code{"eb_zero_null"}; \code{"alra_quantile"} retains the historical detector.}
\item{variance_prior_df}{Equivalent prior sample size for zero-null variance shrinkage.}
\item{factor_target}{Recovery target estimand. Tree-local recovery supports the positive-conditional target only.}
\item{tree_weight}{Relative contribution of available hierarchy proximity to the local kernel, in [0,1].}
\item{tree_tau}{Positive hierarchy-distance decay scale. \code{NULL} estimates a robust scale inside each hard stratum.}
\item{local_k}{Neighbor rank defining the adaptive embedding bandwidth.}
\item{candidate_k}{Embedding-neighbor candidate-pool size. Final-membership cells are also retained as a high-weight core.}
\item{min_effective_donors}{Minimum Kish effective positive-donor size required for tree-local recovery.}
\item{local_info_kappa}{Non-negative information shrinkage strength applied to local residual-factor contribution.}
}
\details{The final SuperCell membership is a cut through a retained walktrap hierarchy, not an assumption that all cells in the cut are exchangeable. In \code{"tree_local_factor"}, hard biological strata are absolute borrowing boundaries. Within a stratum, a query cell can borrow from nearby sibling memberships, but contribution decays according to both hierarchy proximity and adaptive embedding distance. The positive baseline is query-specific, and non-target factors predict only residual expression beyond that baseline using locally weighted ridge regression. Weak local support shrinks residual correction toward the query-specific local mean. The deterministic matrix contains conditional means; DV-sensitive analysis should propagate stored predictive uncertainty or repeated Gamma moment-matched draws.}
\value{A \code{DropoutKillerResult} containing expression, mask, event table, memberships, membership hierarchy when internally built, detection details, settings, uncertainty status, and optional sparse predictive variance. Tree-local results also contain \code{local_geometry}; event tables include local positive mean/variance/prevalence, effective donors, weighted tree/embedding distance diagnostics, factor predictability, shrinkage, and prediction standard deviation.}
