\name{weighted_neighbor_prediction}
\alias{weighted_neighbor_prediction}
\alias{prepare_gene_prior}
\alias{combine_gene_prior}
\alias{gene_prior_from_edges}
\alias{gene_prior_prediction}
\alias{recover_dropout_expression}
\title{Selective cell-space and PPI/pathway-prior recovery components}
\description{Predicts masked dropout events using membership-local latent neighbors and optional PPI/pathway priors, then applies selective hybrid recovery. GRN/TF-target priors are not supported.}
\usage{
weighted_neighbor_prediction(x, embedding, membership, mask, k = 30L,
  sigma = NULL, min_positive_neighbors = 1L, positive_only = TRUE,
  return_events = FALSE)

gene_prior_from_edges(edges, genes = NULL, source = "source", target = "target",
  weight = NULL, prior_type = c("ppi", "pathway"), directed = NULL)

prepare_gene_prior(prior, genes, remove_self = TRUE, normalize = TRUE)

combine_gene_prior(genes, ppi = NULL, pathway = NULL, weights = NULL,
  remove_self = TRUE)

gene_prior_prediction(x, prior, membership, mask, standardize = TRUE,
  min_sd = 1e-8, return_events = FALSE)

recover_dropout_expression(x, mask, membership, embedding, alpha = 0.75,
  gene_prior = NULL, neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  prior_standardize = TRUE, cap_quantile = NULL, return_details = FALSE)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{embedding}{Cell-by-dimension latent representation.}
\item{membership}{One membership label per cell.}
\item{mask}{Sparse or dense logical matrix of eligible zero events.}
\item{k, neighbor_k}{Maximum number of membership-local neighbors.}
\item{sigma, neighbor_sigma}{Gaussian bandwidth. If NULL, an adaptive per-cell median distance is used.}
\item{min_positive_neighbors}{Minimum positive-expression donors required.}
\item{positive_only, neighbor_positive_only}{Use only donors with observed positive expression for the target gene; set FALSE to include zero-valued neighbors in the weighted mean.}
\item{return_events, return_details}{Return event-level details in addition to predictions.}
\item{edges}{Data frame of PPI or pathway edges.}
\item{source, target}{Column names for source and target genes. The matrix orientation is target-by-source.}
\item{weight}{Optional edge-weight column name.}
\item{prior_type}{Either \code{"ppi"} or \code{"pathway"}.}
\item{directed}{Whether edges are directed. If NULL, PPI defaults to symmetric and pathway defaults to directed.}
\item{prior}{One gene-by-gene PPI/pathway prior matrix.}
\item{genes}{Expression gene names.}
\item{remove_self}{Remove self-edges.}
\item{normalize}{L1-normalize target rows.}
\item{ppi}{Optional PPI adjacency prior.}
\item{pathway}{Optional pathway adjacency prior.}
\item{weights}{Optional non-negative PPI/pathway fusion weights.}
\item{standardize, prior_standardize}{Perform membership-local gene scaling before prior prediction.}
\item{min_sd}{Minimum gene standard deviation for standardized prediction.}
\item{alpha}{Cell-space weight when both cell and PPI/pathway predictions are available.}
\item{gene_prior}{Prepared PPI/pathway prior for direct recovery.}
\item{cap_quantile}{Optional upper cap based on positive target-gene expression in the membership.}
}
\value{Sparse predictions, a prepared PPI/pathway prior, or a recovered expression matrix depending on the function.}
