\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Builds or accepts biologically restricted cell memberships, detects selected technical-zero candidates, and selectively recovers only those events. Version 0.7 adds a soft-hierarchy barycentric engine that can use final SuperCell membership as a prior rather than an absolute borrowing wall while preserving explicit hard biological strata. Historical masked-factor, tree-local-factor, and neighbor engines remain available. Observed values are never overwritten.}
\usage{
dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "barycentric", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, membership_penalty = 1,
  barycentric_lambda = 1, barycentric_iter = 20L,
  barycentric_tol = 1e-6)

DropoutKiller(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "barycentric", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, membership_penalty = 1,
  barycentric_lambda = 1, barycentric_iter = 20L,
  barycentric_tol = 1e-6)

run_dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "tree_local_factor", "barycentric", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, membership_penalty = 1,
  barycentric_lambda = 1, barycentric_iter = 20L,
  barycentric_tol = 1e-6)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix, typically normalized/log-transformed expression rather than raw counts.}
\item{embedding}{Cell-by-dimension biological representation used for membership construction and local recovery geometry.}
\item{membership}{Optional membership vector or \code{DropoutKillerMembership}. Internally built memberships retain the full walktrap hierarchy.}
\item{group}{Optional hard biological stratum such as major cell type. Different groups never borrow under tree-local or barycentric recovery.}
\item{split_by}{Optional additional hard boundary such as condition or donor.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters. \code{gamma} controls the final tree cut. Under \code{"barycentric"}, the cut is used as a soft prior rather than an absolute borrowing boundary.}
\item{rank, max_rank, rank_z, quantile_prob, threshold, min_cells, min_negative}{Dropout-detection controls.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor controls; positional slots are retained.}
\item{cap_quantile}{Optional historical target-gene cap.}
\item{seed}{Random seed.}
\item{return_score}{Also return the sparse detection-score matrix.}
\item{recovery_method}{Recovery engine. \code{"masked_factor"} remains the default comparator. \code{"tree_local_factor"} uses hierarchy- and embedding-weighted positive baselines inside final memberships plus residual factors. \code{"barycentric"} learns target-independent simplex donor weights across embedding-local cells within hard strata while using membership and hierarchy as soft priors. \code{"neighbor"} is the legacy Gaussian cell-neighbor engine.}
\item{factor_rank, factor_features, factor_ridge}{Non-target coexpression-factor and ridge parameters used by factor engines.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum positive target donors before factor residual learning is attempted.}
\item{detection_method}{Detection engine. Default \code{"eb_zero_null"}; \code{"alra_quantile"} retains the historical detector.}
\item{variance_prior_df}{Equivalent prior sample size for zero-null variance shrinkage.}
\item{factor_target}{Recovery target estimand. Tree-local and barycentric recovery support the positive-conditional target only.}
\item{tree_weight}{Relative contribution of available hierarchy proximity to the local geometry prior, in [0,1].}
\item{tree_tau}{Positive hierarchy-distance decay scale. \code{NULL} estimates a robust scale from available local tree distances.}
\item{local_k}{Neighbor rank defining the adaptive embedding bandwidth.}
\item{candidate_k}{Maximum embedding-neighbor candidate-pool size. Barycentric candidates may cross final membership inside the same hard stratum.}
\item{min_effective_donors}{Minimum Kish effective positive-donor size required for tree-local or barycentric recovery.}
\item{local_info_kappa}{Non-negative information shrinkage strength used by the tree-local residual-factor engine.}
\item{membership_penalty}{Non-negative additive cost for a barycentric candidate outside the query cell's final membership. Zero removes the membership preference; larger values increasingly favor same-membership donors without making the boundary absolute.}
\item{barycentric_lambda}{Positive penalty shrinking barycentric simplex weights toward the hierarchy/embedding/membership donor prior.}
\item{barycentric_iter}{Maximum projected-gradient iterations for each query-cell barycentric weight vector.}
\item{barycentric_tol}{Positive convergence tolerance for barycentric projected-gradient updates.}
}
\details{In \code{"barycentric"}, explicit hard biological strata remain absolute borrowing boundaries. Within a hard stratum, embedding k-nearest neighbors define the local candidate set. Embedding distance, retained SuperCell hierarchy when available, and \code{membership_penalty} define a donor prior. A convex simplex-constrained objective then learns one query-cell weight vector that reconstructs the query in biological-embedding space while remaining close to that prior. The target gene does not participate in weight learning. For each dropout gene, the cell-state weights are renormalized over reliable positive target-gene donors and used to compute a positive-conditional local mean. Kish effective donor size and weighted biological variance are retained for uncertainty diagnostics. Existing factor and neighbor engines are unchanged so recovery families can be compared on identical oracle masks.}
\value{A \code{DropoutKillerResult} containing expression, mask, event table, memberships, membership hierarchy when internally built, detection details, settings, uncertainty status, and optional sparse predictive variance. Barycentric results contain query-cell sparse donor geometry and event-level local positive mean/variance/prevalence, effective donors, weighted tree/embedding distance diagnostics, state-reconstruction improvement, and prediction standard deviation.}
