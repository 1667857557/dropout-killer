\name{dropout_killer}
\alias{dropout_killer}
\alias{DropoutKiller}
\alias{run_dropout_killer}
\title{Run the complete DropoutKiller workflow}
\description{Builds or accepts cell memberships, detects high-confidence technical zero candidates with local ALRA-derived thresholding, and selectively recovers them using membership-local Gaussian neighbor weights estimated from latent-space distances.}
\usage{
dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 20, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE)

DropoutKiller(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 20, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE)

run_dropout_killer(x, embedding, membership = NULL, group = NULL,
  split_by = NULL, gamma = 20, k_knn = 5L, approximate = "auto",
  approx_n = 20000L, rank = "auto", max_rank = 20L, rank_z = 6,
  quantile_prob = 0.001, threshold = 0.95, min_cells = 8L,
  min_negative = 3L, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, seed = 12345L, return_score = FALSE)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix, preferably normalized/log-transformed expression rather than raw counts.}
\item{embedding}{Cell-by-dimension latent representation used for membership and neighbor distances.}
\item{membership}{Optional precomputed memberships.}
\item{group}{Optional major-cell-type or other hard biological stratum.}
\item{split_by}{Optional additional hard boundary such as condition or donor.}
\item{gamma, k_knn, approximate, approx_n}{Membership-construction parameters.}
\item{rank, max_rank, rank_z, quantile_prob, threshold, min_cells, min_negative}{Local low-rank detection parameters.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors}{Neighbor-weight recovery parameters. With NULL bandwidth, the per-cell bandwidth is estimated from neighbor distances.}
\item{neighbor_positive_only}{If TRUE, renormalize Gaussian weights over neighbors with observed positive target-gene expression; FALSE uses the all-neighbor weighted mean.}
\item{cap_quantile}{Optional membership-local positive-expression cap.}
\item{seed}{Random seed.}
\item{return_score}{Also return the sparse score matrix.}
}
\value{A \code{DropoutKillerResult} containing the recovered expression matrix, sparse mask, event table, memberships, detection details, and settings.}
