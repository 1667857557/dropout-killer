\name{weighted_neighbor_prediction}
\alias{weighted_neighbor_prediction}
\alias{recover_dropout_expression}
\title{Membership-local weighted-neighbor recovery components}
\description{Predicts masked dropout events using membership-constrained latent-space neighbors and applies selective zero-only recovery. Donor weights are estimated from Gaussian functions of latent-space distance; no external biological prior is used.}
\usage{
weighted_neighbor_prediction(x, embedding, membership, mask, k = 30L,
  sigma = NULL, min_positive_neighbors = 1L, positive_only = TRUE,
  return_events = FALSE)

recover_dropout_expression(x, mask, membership, embedding,
  neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, return_details = FALSE)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{embedding}{Cell-by-dimension latent representation.}
\item{membership}{One membership label per cell.}
\item{mask}{Sparse or dense logical matrix of eligible zero events.}
\item{k, neighbor_k}{Maximum number of membership-local neighbors.}
\item{sigma, neighbor_sigma}{Gaussian bandwidth. If NULL, an adaptive per-cell median positive neighbor distance is used.}
\item{min_positive_neighbors}{Minimum positive-expression donors required.}
\item{positive_only, neighbor_positive_only}{Use only donors with observed positive expression for the target gene; set FALSE to include zero-valued neighbors in the weighted mean.}
\item{return_events, return_details}{Return event-level details in addition to predictions.}
\item{cap_quantile}{Optional upper cap based on positive target-gene expression in the membership.}
}
\value{A sparse prediction matrix, an event table when requested, or a selectively recovered expression matrix.}
