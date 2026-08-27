\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Applies ALRA-style library-size and log normalization by default, detects selected technical-zero candidates, and selectively recovers only those events. The production detector uses native-ALRA-style randomized low-rank reconstruction and the gene-wise empirical 0.1 percent quantile zero gate globally within each supplied major cell class. SuperCell memberships and retained hierarchy remain recovery structures rather than zero-detection blocks.}
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
  detection_method = c("alra_global_by_group", "eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, alra_K = 100L,
  alra_noise_start = 80L, alra_choose_q = 2L, alra_svd_q = 10L)

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
  detection_method = c("alra_global_by_group", "eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, alra_K = 100L,
  alra_noise_start = 80L, alra_choose_q = 2L, alra_svd_q = 10L)

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
  detection_method = c("alra_global_by_group", "eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, alra_K = 100L,
  alra_noise_start = 80L, alra_choose_q = 2L, alra_svd_q = 10L)
}
\arguments{
\item{x}{Non-negative gene-by-cell raw count/expression matrix. By default it is internally library-size normalized and log transformed using the ALRA convention. Set \code{normalize=FALSE} only when supplying an already normalized working matrix.}
\item{embedding}{Cell-by-dimension biological representation used for membership construction and local recovery geometry.}
\item{membership}{Optional membership vector or \code{DropoutKillerMembership}. Membership controls recovery borrowing; it does not fragment the default global-ALRA detection blocks.}
\item{group}{Optional major cell-class vector. Under the default detector each group is one independent global ALRA block. If \code{NULL}, all cells form one ALRA block.}
\item{split_by}{Optional additional hard recovery boundary such as condition or donor. It does not split default ALRA detection blocks.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters; \code{gamma} controls the final tree cut but the retained hierarchy remains available to recovery.}
\item{rank}{ALRA rank. \code{"auto"} uses the original singular-value-spacing heuristic independently in each cell-class block; a numeric value fixes rank subject to block dimensions.}
\item{max_rank}{Maximum rank used only by the historical membership-local detectors.}
\item{rank_z}{Spacing threshold in standard deviations for automatic rank selection. Default 6, matching ALRA.}
\item{quantile_prob}{Per-gene low-rank quantile used by ALRA zero detection. Default 0.001.}
\item{threshold}{Selection threshold for historical confidence-based detectors. It is not applied to \code{"alra_global_by_group"}, whose native ALRA quantile rule is already binary.}
\item{min_cells}{Minimum number of cells required for a detection block.}
\item{min_negative}{Minimum negative low-rank support used only by historical membership-local detectors.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor recovery controls.}
\item{cap_quantile}{Optional historical target-gene cap.}
\item{seed}{Random seed.}
\item{return_score}{Also return a sparse detection-score matrix. For global ALRA this is a binary call matrix.}
\item{recovery_method}{Recovery engine. \code{"masked_factor"} is the default comparator. \code{"tree_local_factor"} uses retained hierarchy plus embedding-local positive baselines and locally weighted residual factors. \code{"neighbor"} is the legacy cell-neighbor engine.}
\item{factor_rank, factor_features, factor_ridge}{Non-target coexpression-factor and ridge parameters.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum positive target donors before factor residual learning is attempted.}
\item{detection_method}{Detection engine. Default \code{"alra_global_by_group"} runs native-ALRA-style global detection independently within each \code{group}. \code{"eb_zero_null"} and \code{"alra_quantile"} retain the historical membership-local comparators.}
\item{variance_prior_df}{Equivalent prior sample size for the historical EB zero-null variance shrinkage.}
\item{factor_target}{Recovery target estimand. Tree-local recovery supports the positive-conditional target only.}
\item{tree_weight}{Relative contribution of available hierarchy proximity to the local kernel, in [0,1].}
\item{tree_tau}{Positive hierarchy-distance decay scale. \code{NULL} estimates a robust scale inside each hard recovery stratum.}
\item{local_k}{Neighbor rank defining the adaptive embedding bandwidth.}
\item{candidate_k}{Embedding-neighbor candidate-pool size within the final membership.}
\item{min_effective_donors}{Minimum Kish effective positive-donor size required for tree-local recovery.}
\item{local_info_kappa}{Non-negative information shrinkage strength applied to local residual-factor contribution.}
\item{normalize}{Logical. If \code{TRUE} (default), apply ALRA-equivalent library-size normalization followed by \code{log1p} once before detection and recovery.}
\item{normalization_scale_factor}{Positive library-size target used before \code{log1p}; default \code{1e4}, matching ALRA.}
\item{alra_K}{Maximum number of singular values inspected by automatic global-ALRA rank selection. Default 100, matching original ALRA when block dimensions permit.}
\item{alra_noise_start}{Original ALRA noise-tail start index. Default 80; automatically reduced for blocks with fewer than 100 available singular values while retaining at least five tail spacings when possible.}
\item{alra_choose_q}{Randomized-SVD power iterations for rank selection. Default 2, matching original ALRA.}
\item{alra_svd_q}{Randomized-SVD power iterations for final rank-k reconstruction. Default 10, matching original ALRA.}
}
\details{For genes in rows and cells in columns, default preprocessing is \deqn{X'_{gc}=\log(1+10^4 X_{gc}/\sum_h X_{hc}).} For each major cell class, the default detector transposes the class block to the original ALRA cell-by-gene orientation, chooses rank from singular-value spacings when \code{rank="auto"}, reconstructs the class matrix, and defines \deqn{\tau_g=|Q_{0.001}(\hat X_{\cdot g})|.} An observed zero is selected when its low-rank reconstruction is strictly greater than \eqn{\tau_g}. Final SuperCell membership, retained hierarchy, \code{split_by}, and embedding-local geometry are then used only by recovery. This separation prevents small SuperCell memberships from imposing the finite-sample detection bottleneck observed in the zero-detection benchmark.}
\value{A \code{DropoutKillerResult} containing normalized recovered expression, mask, event table, memberships, membership hierarchy when internally built, detection details, settings, uncertainty status, and optional sparse predictive variance. Global-ALRA events additionally contain \code{detection_block} and \code{alra_margin}.}
