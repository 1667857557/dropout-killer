\name{local_alra_detect}
\alias{local_alra_detect}
\alias{local_alra_score}
\alias{select_dropout_mask}
\title{Membership-local ALRA-inspired dropout detection}
\description{Computes low-rank reconstructions independently within memberships, applies a gene-specific negative-tail ALRA gate, and assigns a null-tail confidence score only to observed zero events.}
\usage{
local_alra_detect(x, membership, rank = "auto", max_rank = 20L,
  rank_z = 6, quantile_prob = 0.001, min_cells = 8L,
  min_negative = 3L, seed = 1L)

local_alra_score(x, membership, rank = "auto", max_rank = 20L,
  rank_z = 6, quantile_prob = 0.001, min_cells = 8L,
  min_negative = 3L, seed = 1L)

select_dropout_mask(score, threshold = 0.95, x = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell.}
\item{rank}{Positive integer or \code{"auto"}.}
\item{max_rank}{Maximum rank considered in automatic selection.}
\item{rank_z}{Singular-spacing signal threshold in tail standard deviations.}
\item{quantile_prob}{Lower quantile used for the ALRA-derived adaptive threshold.}
\item{min_cells}{Minimum membership size for local decomposition.}
\item{min_negative}{Minimum negative reconstructed values for parametric confidence estimation.}
\item{seed}{Random seed used by randomized SVD when required.}
\item{score}{Detection object or score matrix.}
\item{threshold}{High-confidence cutoff.}
}
\details{The confidence is evidence against a symmetric biological-zero reconstruction-error null and is not a calibrated posterior probability of dropout.}
\value{A detection object, sparse confidence matrix, or sparse logical mask depending on the function.}
