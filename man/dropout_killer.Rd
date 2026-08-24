\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Builds or accepts biologically restricted cell memberships, detects observed zeros whose membership-local low-rank reconstruction is inconsistent with a finite-sample symmetric biological-zero null, and selectively recovers only selected events. The default detector uses empirical-Bayes negative-tail variance shrinkage with gene-wise Benjamini-Hochberg adjustment. The default recovery engine estimates positive conditional target magnitude from membership-local non-target coexpression factors with ridge regression, analytic leave-one-out shrinkage, and predictive uncertainty. Observed values are never overwritten.}
\usage{
dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"))

DropoutKiller(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"))

run_dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10,
  factor_target = c("positive", "all_observed"))
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix, typically normalized/log-transformed expression rather than raw counts.}
\item{embedding}{Cell-by-dimension representation used for membership construction and by the optional neighbor recovery engine.}
\item{membership}{Optional precomputed memberships. Internally built memberships are canonically relabeled before detection/recovery without changing the partition.}
\item{group}{Optional hard biological stratum such as major cell type.}
\item{split_by}{Optional additional hard boundary such as condition or donor.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters; \code{gamma} defaults to 150.}
\item{rank, max_rank, rank_z}{Membership-local low-rank reconstruction and automatic-rank parameters.}
\item{quantile_prob}{Historical ALRA lower-tail quantile, used only with \code{detection_method="alra_quantile"}.}
\item{threshold}{Detection evidence cutoff. Under the default detector, confidence equals \code{1-q_value}; 0.95 therefore selects approximately gene-wise BH q <= 0.05 under the working zero-null.}
\item{min_cells}{Minimum membership size for local low-rank detection.}
\item{min_negative}{Negative-reconstruction support parameter for zero-null variance estimation.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor recovery controls. Their existing positional slots are retained; they are used only with \code{recovery_method="neighbor"}.}
\item{cap_quantile}{Optional target-gene membership-local upper cap.}
\item{seed}{Random seed.}
\item{return_score}{Also return the sparse detection-score matrix.}
\item{recovery_method}{Recovery engine; default \code{"masked_factor"}.}
\item{factor_rank, factor_features, factor_ridge}{Coexpression-factor and ridge settings. Factor learning uses a truncated SVD of standardized non-target genes.}
\item{min_feature_observed}{Minimum observations for a non-target factor-learning feature.}
\item{min_target_observed}{Minimum eligible target-gene donors for factor regression.}
\item{detection_method}{Detection engine. Default \code{"eb_zero_null"} uses finite-sample negative-tail variance shrinkage and gene-wise BH testing. \code{"alra_quantile"} retains the historical empirical extreme-quantile gate.}
\item{variance_prior_df}{Non-negative equivalent prior sample size used to shrink gene-specific negative-tail zero-null variance toward the robust membership-level center.}
\item{factor_target}{Masked-factor target estimand. Default \code{"positive"} learns target magnitude from reliable positive donors after a zero has already been selected as technical dropout. \code{"all_observed"} reproduces the previous zero-inflated target.}
}
\details{The default detector's \code{confidence} is an evidence score \code{1-q_value}, not a calibrated Bayesian posterior probability of technical dropout. Gene-wise BH control is conditional on the approximate symmetric Gaussian low-rank zero-null and is not claimed to be exact global FDR control across genes or memberships. Positive-conditional recovery must only be applied after selective dropout detection; it is not appropriate for indiscriminate imputation of natural zeros. The deterministic completed matrix contains conditional means and does not by itself preserve full differential variability; use event predictive variance or repeated completed draws for DV-sensitive analyses.}
\value{A \code{DropoutKillerResult} containing the deterministic recovered expression matrix, mask, event table, memberships, detection details, settings, and \code{uncertainty_available}. Event tables from the default detector include zero-null scale, z score, p value, q value, and confidence. For masked-factor results, \code{predictive_variance} is a sparse matrix when complete event-level predictive uncertainty is available. For neighbor results it is \code{NULL}.}
