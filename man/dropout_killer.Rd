\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{
Detects selected technical-zero candidates and selectively recovers only those
coordinates. The default detector is a source-faithful port of the scGACL
Gamma-Normal mixture detector operating on raw counts within supplied cell
subpopulations. Original all-cell ALRA remains an explicit comparator. The
production recovery engine is the target-safe P1 stabilized state.
}
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
  detection_method = c("scgacl_gamma_normal", "alra_global",
  "alra_global_by_group", "eb_zero_null", "alra_quantile"),
  scgacl_dropout_threshold = 0.5, variance_prior_df = 10,
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
  detection_method = c("scgacl_gamma_normal", "alra_global",
  "alra_global_by_group", "eb_zero_null", "alra_quantile"),
  scgacl_dropout_threshold = 0.5, variance_prior_df = 10,
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
  detection_method = c("scgacl_gamma_normal", "alra_global",
  "alra_global_by_group", "eb_zero_null", "alra_quantile"),
  scgacl_dropout_threshold = 0.5, variance_prior_df = 10,
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
\item{x}{Non-negative gene-by-cell matrix. For the default scGACL detector this must be the raw count matrix.}
\item{embedding}{Cell-by-dimension biological representation used for membership construction and local recovery geometry.}
\item{membership}{Optional membership vector or \code{DropoutKillerMembership}. Membership controls recovery borrowing and does not fragment the original global ALRA comparator.}
\item{group}{Cell-subpopulation labels. Required by the default source-faithful scGACL path because the released scGACL configuration uses known cell-type labels. Also acts as a hard recovery boundary.}
\item{split_by}{Optional additional hard recovery boundary such as condition or donor.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters.}
\item{rank}{ALRA rank. \code{"auto"} uses the original singular-value-spacing heuristic; a positive numeric value fixes rank.}
\item{max_rank}{Maximum rank used only by historical membership-local detectors.}
\item{rank_z}{Spacing threshold in standard deviations for original ALRA automatic rank selection. Default 6.}
\item{quantile_prob}{Per-gene low-rank quantile used by ALRA zero detection. Default 0.001.}
\item{threshold}{Selection threshold for historical confidence-based local detectors. It is not an additional gate for scGACL or original ALRA.}
\item{min_cells}{Minimum membership size used by historical local detectors.}
\item{min_negative}{Minimum negative low-rank support used only by historical membership-local detectors.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor recovery controls.}
\item{cap_quantile}{Optional historical target-gene cap.}
\item{seed}{Random seed. For ALRA, one stream is used across rank selection and final randomized SVD as in the released workflow.}
\item{return_score}{Also return a sparse selected-event score matrix. scGACL scores are Gamma-component dropout posteriors; ALRA scores are binary calls.}
\item{recovery_method}{Recovery engine. \code{"p1_stabilized_state"} is the production default.}
\item{factor_rank, factor_features, factor_ridge}{Non-target factor-state and ridge parameters.}
\item{min_feature_observed}{Minimum observations for a non-target factor feature.}
\item{min_target_observed}{Minimum positive target donors before factor recovery is attempted.}
\item{detection_method}{Detection engine. Default \code{"scgacl_gamma_normal"} uses the released scGACL Gamma-Normal mixture detector. \code{"alra_global"} uses original all-cell ALRA. \code{"alra_global_by_group"} is a deprecated compatibility alias that warns and maps to \code{"alra_global"}. \code{"eb_zero_null"} and \code{"alra_quantile"} retain historical membership-local comparators.}
\item{scgacl_dropout_threshold}{Gamma-component posterior threshold for the scGACL detector. Default 0.5, matching the released configuration.}
\item{variance_prior_df}{Equivalent prior sample size for historical EB zero-null variance shrinkage.}
\item{factor_target}{Recovery target estimand.}
\item{tree_weight}{Relative hierarchy contribution to the recovery kernel, in [0,1].}
\item{tree_tau}{Positive hierarchy-distance decay scale.}
\item{local_k}{Neighbor rank defining adaptive recovery bandwidth.}
\item{candidate_k}{Embedding-neighbor candidate-pool size inside final membership.}
\item{min_effective_donors}{Minimum Kish effective positive-donor size for tree-local recovery.}
\item{local_info_kappa}{Non-negative information shrinkage strength for local residual-factor contribution.}
\item{normalize}{Logical. If true, construct the recovery working matrix by ALRA library-size normalization followed by \code{log1p}. This does not alter the raw-count input used by scGACL detection.}
\item{normalization_scale_factor}{Positive recovery library-size target; default \code{1e4}.}
\item{alra_K}{Number of singular values inspected by original ALRA automatic rank selection. Default 100. The original source validity check is retained rather than adapting K to small blocks.}
\item{alra_noise_start}{Original ALRA noise-tail start index. Default 80.}
\item{alra_choose_q}{Randomized-SVD power iterations for ALRA rank selection. Default 2.}
\item{alra_svd_q}{Randomized-SVD power iterations for final ALRA reconstruction. Default 10.}
\item{factor_crossfit_folds}{Number of deterministic target-gene folds used for leakage-safe factor-state construction.}
\item{factor_crossfit_seed}{Integer controlling deterministic target-gene fold assignment.}
\item{support_adaptive_rank}{Logical; cap target regression rank using positive-donor support.}
\item{bias_kappa}{Non-negative LOO bias-shrinkage constant.}
\item{predictor_smoothing}{Weight in \eqn{[0,1]} for the single hierarchy/embedding smoothing step applied only to the target-safe predictor state.}
}
\details{
For the default scGACL detector, raw counts are transformed as
\deqn{x_{gc}=\log(1.01+C_{gc}).}
Within each supplied cell subpopulation and gene, the released model fits
\deqn{f(x)=\lambda\,Gamma(x;\alpha,\beta)+(1-\lambda)N(x;\mu,\sigma),}
where \eqn{\beta} is a rate. The dropout score is the Gamma responsibility
\deqn{d(x)=\frac{\lambda Gamma(x;\alpha,\beta)}{\lambda Gamma(x;\alpha,\beta)+(1-\lambda)N(x;\mu,\sigma)}.}
Only raw-zero coordinates are eligible and the default released threshold is
0.5.

Detection and recovery scales are independent. With \code{normalize=TRUE},
recovery uses
\deqn{X'_{gc}=\log(1+10^4 X_{gc}/\sum_h X_{hc}).}

For \code{detection_method="alra_global"}, the complete all-cell recovery-scale
matrix is transposed to the original ALRA cell-by-gene orientation. Automatic
rank uses the released singular-value-spacing rule with \code{K=100},
\code{noise_start=80}, threshold 6 and \code{q=2}. The final randomized SVD uses
\code{q=10} and defines
\deqn{\tau_g=|Q_{0.001}(\hat X_{\cdot g})|.}
An observed zero is selected when \eqn{\hat X_{cg}>\tau_g}. No cell-type or
membership-specific ALRA SVD is fitted.
}
\value{
A \code{DropoutKillerResult} containing the recovered expression working
matrix, sparse mask, event table, memberships, detection details, settings,
uncertainty status, and optional predictive variance. scGACL events retain the
mixture posterior and fitted Gamma/Normal parameters; ALRA events retain low-rank
value, threshold and margin.
}
