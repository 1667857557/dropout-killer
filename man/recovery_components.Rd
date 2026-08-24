\name{masked_factor_prediction}
\alias{masked_factor_prediction}
\alias{weighted_neighbor_prediction}
\alias{recover_dropout_expression}
\alias{sample_dropout_expression}
\title{Selective dropout recovery components}
\description{Predicts supplied dropout-mask events with membership-local masked coexpression factors and ridge regression. Target genes are excluded from factor-state learning and analytic leave-one-out predictions choose shrinkage toward the membership mean. The legacy Gaussian neighbor estimator remains available explicitly. Predictive uncertainty is retained for differential-variability-aware downstream analysis.}
\usage{
masked_factor_prediction(x, membership, mask, factor_rank = 5L,
  factor_features = 2000L, factor_max_iter = 8L,
  factor_tol = 1e-04, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  cap_quantile = NULL, return_events = FALSE)

weighted_neighbor_prediction(x, embedding, membership, mask, k = 30L,
  sigma = NULL, min_positive_neighbors = 1L, positive_only = TRUE,
  return_events = FALSE)

recover_dropout_expression(x, mask, membership, embedding = NULL,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L,
  factor_max_iter = 8L, factor_tol = 1e-04, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, return_details = FALSE)

sample_dropout_expression(result, n = 1L, seed = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell.}
\item{mask}{Sparse or dense logical matrix containing only zero-valued dropout events.}
\item{embedding}{Optional cell-by-dimension latent representation. Required only for \code{recovery_method="neighbor"}.}
\item{recovery_method}{Recovery engine. The default \code{"masked_factor"} treats mask entries as missing. \code{"neighbor"} selects the legacy Gaussian cell-neighbor estimator.}
\item{factor_rank}{Maximum membership-local coexpression-factor rank.}
\item{factor_features}{Maximum number of high-variance non-target genes used to learn membership cell factors.}
\item{factor_max_iter}{Maximum masked low-rank update iterations.}
\item{factor_tol}{Relative convergence tolerance for masked factor updates.}
\item{factor_ridge}{Non-negative ridge penalty for target-gene factor coefficients; the intercept is unpenalized.}
\item{min_feature_observed}{Minimum unmasked observations required for a non-target gene to participate in factor learning.}
\item{min_target_observed}{Minimum unmasked target-gene observations for factor regression; otherwise the membership mean is used.}
\item{cap_quantile}{Optional upper cap based on the target gene's unmasked membership expression.}
\item{return_events, return_details}{Return event-level diagnostics in addition to predictions/recovered expression.}
\item{k, neighbor_k}{Maximum membership-local neighbors for the legacy neighbor engine.}
\item{sigma, neighbor_sigma}{Gaussian bandwidth for the legacy neighbor engine.}
\item{min_positive_neighbors}{Minimum positive-expression donors required by the legacy neighbor engine.}
\item{positive_only, neighbor_positive_only}{Whether the legacy neighbor estimator excludes zero-valued donors.}
\item{result}{A \code{DropoutKillerResult} containing event-level predictive standard deviations.}
\item{n}{Number of uncertainty-aware completed matrices to draw.}
\item{seed}{Optional random seed for completed-matrix draws.}
}
\value{\code{masked_factor_prediction} and \code{weighted_neighbor_prediction} return sparse predictions unless event tables are requested. \code{recover_dropout_expression} returns a selectively recovered matrix or a list with event diagnostics. \code{sample_dropout_expression} returns one completed matrix when \code{n=1} and a list otherwise.}
