\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Builds or accepts cell memberships, detects high-confidence technical zero candidates with membership-local ALRA-inspired thresholding, and selectively recovers them. The default recovery engine excludes target genes from factor-state learning and uses membership-local coexpression factors with target-gene ridge regression, analytic leave-one-out shrinkage, and predictive uncertainty.}
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
  min_feature_observed = 20L, min_target_observed = 20L)

DropoutKiller(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L)

run_dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 150, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix, typically normalized/log-transformed expression rather than raw counts.}
\item{embedding}{Cell-by-dimension representation used for membership construction and by the optional neighbor recovery engine.}
\item{membership}{Optional precomputed memberships.}
\item{group}{Optional hard biological stratum such as major cell type.}
\item{split_by}{Optional additional hard boundary such as condition or donor.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters; \code{gamma} defaults to 150.}
\item{rank, max_rank, rank_z, quantile_prob, threshold, min_cells, min_negative}{Membership-local dropout-detection parameters.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Historical Gaussian-neighbor recovery controls. Their pre-0.4 positional slots are retained; they are used only with \code{recovery_method="neighbor"}.}
\item{cap_quantile}{Optional target-gene membership-local upper cap.}
\item{seed}{Random seed.}
\item{return_score}{Also return the sparse detection-score matrix.}
\item{recovery_method}{Recovery engine; default \code{"masked_factor"}. New recovery-engine controls are appended after the historical positional API.}
\item{factor_rank, factor_features, factor_ridge}{Coexpression-factor and ridge settings. Factor learning uses a truncated SVD of standardized non-target genes.}
\item{min_feature_observed}{Minimum unmasked observations for a non-target factor-learning feature.}
\item{min_target_observed}{Minimum unmasked target-gene observations for factor regression.}
}
\value{A \code{DropoutKillerResult} containing the deterministic recovered expression matrix, mask, event table, memberships, detection details, settings, and \code{uncertainty_available}. For masked-factor results, \code{predictive_variance} is a sparse matrix. For neighbor results it is \code{NULL}, because uncertainty is unavailable rather than zero.}
