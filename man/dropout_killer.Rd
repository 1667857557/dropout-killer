\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Applies ALRA-style library-size and log normalization by default, detects selected technical-zero candidates, and selectively recovers only those events. The production detector uses native-ALRA-style randomized low-rank reconstruction and the gene-wise empirical 0.1 percent quantile zero gate globally within each supplied major cell class. The production imputation engine is the target-safe P1 stabilized state.}
\usage{
dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("p1_stabilized_state", "masked_factor",
  "tree_local_factor", "neighbor"), factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 2,
  min_feature_observed = 20L, min_target_observed = 8L,
  detection_method = c("alra_global_by_group", "eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, alra_K = 100L,
  alra_noise_start = 80L, alra_choose_q = 2L, alra_svd_q = 10L,
  factor_crossfit_folds = 5L, factor_crossfit_seed = 1L,
  support_adaptive_rank = TRUE, bias_kappa = 10,
  predictor_smoothing = 0.25)

DropoutKiller(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("p1_stabilized_state", "masked_factor",
  "tree_local_factor", "neighbor"), factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 2,
  min_feature_observed = 20L, min_target_observed = 8L,
  detection_method = c("alra_global_by_group", "eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, alra_K = 100L,
  alra_noise_start = 80L, alra_choose_q = 2L, alra_svd_q = 10L,
  factor_crossfit_folds = 5L, factor_crossfit_seed = 1L,
  support_adaptive_rank = TRUE, bias_kappa = 10,
  predictor_smoothing = 0.25)

run_dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("p1_stabilized_state", "masked_factor",
  "tree_local_factor", "neighbor"), factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 2,
  min_feature_observed = 20L, min_target_observed = 8L,
  detection_method = c("alra_global_by_group", "eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"),
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, normalize = TRUE,
  normalization_scale_factor = 1e4, alra_K = 100L,
  alra_noise_start = 80L, alra_choose_q = 2L, alra_svd_q = 10L,
  factor_crossfit_folds = 5L, factor_crossfit_seed = 1L,
  support_adaptive_rank = TRUE, bias_kappa = 10,
  predictor_smoothing = 0.25)
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
\item{recovery_method}{Recovery engine. \code{"p1_stabilized_state"} is the production default: target-gene folds are excluded from predictor construction, the predictor state is smoothed once over the hierarchy/embedding geometry, and P1 positive-donor ridge predictions are LOO bias-calibrated. The other values retain the comparison engines.}
\item{factor_rank, factor_features, factor_ridge}{Non-target coexpression-factor and ridge parameters. The production P1 defaults are rank 5 and ridge 2.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum positive target donors before factor recovery is attempted. The production P1 default is 8.}
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
\item{factor_crossfit_folds}{Number of deterministic target-gene folds used for leakage-safe factor-state construction. The production P1 default is 5.}
\item{factor_crossfit_seed}{Integer controlling deterministic target-gene fold assignment.}
\item{support_adaptive_rank}{Logical. If true, cap target regression rank using positive-donor support. Enabled by default for production P1 recovery.}
\item{bias_kappa}{Non-negative LOO bias-shrinkage constant. The production P1 default is 10; \code{Inf} disables calibration.}
\item{predictor_smoothing}{Weight in \eqn{[0,1]} for the single row-stochastic hierarchy/embedding smoothing step applied only to the target-safe predictor state. Default 0.25.}
}
\details{For genes in rows and cells in columns, default preprocessing is \deqn{X'_{gc}=\log(1+10^4 X_{gc}/\sum_h X_{hc}).} For each major cell class, the default detector transposes the class block to the original ALRA cell-by-gene orientation, chooses rank from singular-value spacings when \code{rank="auto"}, reconstructs the class matrix, and defines \deqn{\tau_g=|Q_{0.001}(\hat X_{\cdot g})|.} An observed zero is selected when its low-rank reconstruction is strictly greater than \eqn{\tau_g}. For P1 stabilized-state recovery, each target fold is removed before predictor-state construction, the standardized non-target state receives one row-stochastic geometry smoothing step, and the resulting factors feed positive-donor ridge prediction with analytic LOO shrinkage and bias calibration. Observed non-target coordinates remain exactly unchanged.}
\value{A \code{DropoutKillerResult} containing normalized recovered expression, mask, event table, memberships, membership hierarchy when internally built, detection details, settings, uncertainty status, and optional sparse predictive variance. Global-ALRA events additionally contain \code{detection_block} and \code{alra_margin}.}
