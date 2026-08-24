\name{masked_factor_prediction}
\alias{masked_factor_prediction}
\alias{weighted_neighbor_prediction}
\alias{recover_dropout_expression}
\alias{sample_dropout_expression}
\title{Selective dropout recovery components}
\description{Predicts supplied dropout-mask events with membership-local coexpression factors and ridge regression. Target genes are excluded from factor-state learning. By default, after a coordinate has already been classified as technical dropout, target magnitude is learned only from reliable positive donors. Analytic leave-one-out predictions select shrinkage toward the corresponding membership mean and provide approximate predictive uncertainty. The legacy Gaussian neighbor estimator remains available explicitly.}
\usage{
masked_factor_prediction(x, membership, mask, factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  cap_quantile = NULL, return_events = FALSE,
  target_mode = c("positive", "all_observed"))

weighted_neighbor_prediction(x, embedding, membership, mask, k = 30L,
  sigma = NULL, min_positive_neighbors = 1L, positive_only = TRUE,
  return_events = FALSE)

recover_dropout_expression(x, mask, membership, embedding = NULL,
  neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, return_details = FALSE,
  recovery_method = c("masked_factor", "neighbor"),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 1,
  min_feature_observed = 20L, min_target_observed = 20L,
  factor_target = c("positive", "all_observed"))

sample_dropout_expression(result, n = 1L, seed = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell.}
\item{mask}{Sparse or dense logical matrix containing only zero-valued dropout events.}
\item{embedding}{Optional cell-by-dimension latent representation. Required only for \code{recovery_method="neighbor"}.}
\item{neighbor_k}{Maximum membership-local neighbors for the legacy neighbor engine. This retains its historical positional slot.}
\item{neighbor_sigma}{Gaussian bandwidth for the legacy neighbor engine. This retains its historical positional slot.}
\item{min_positive_neighbors}{Minimum positive-expression donors required by the legacy neighbor engine.}
\item{neighbor_positive_only}{Whether the legacy neighbor estimator excludes zero-valued donors.}
\item{cap_quantile}{Optional upper cap based on the target gene's eligible membership-local donor expression.}
\item{return_details}{Return event-level diagnostics in addition to recovered expression.}
\item{recovery_method}{Recovery engine. \code{"masked_factor"} is the default; \code{"neighbor"} selects the legacy Gaussian cell-neighbor estimator.}
\item{factor_rank}{Maximum membership-local coexpression-factor rank.}
\item{factor_features}{Maximum number of high-variance non-target genes used to learn membership cell factors.}
\item{factor_ridge}{Non-negative ridge penalty for target-gene factor coefficients; the intercept is unpenalized.}
\item{min_feature_observed}{Minimum observations required for a non-target gene to participate in factor learning.}
\item{min_target_observed}{Minimum eligible target-gene donor observations for factor regression; otherwise the engine falls back to the matching membership mean when it can estimate one.}
\item{target_mode}{Target estimand for \code{masked_factor_prediction}. \code{"positive"} estimates expression magnitude using only reliable positive donors after a mask has already classified the query zero as technical dropout. \code{"all_observed"} reproduces the previous zero-inflated target model.}
\item{factor_target}{End-to-end/direct-recovery analogue of \code{target_mode}; default \code{"positive"}.}
\item{return_events}{Return event-level diagnostics from a prediction function.}
\item{k}{Maximum membership-local neighbors for \code{weighted_neighbor_prediction}.}
\item{sigma}{Gaussian bandwidth for \code{weighted_neighbor_prediction}.}
\item{positive_only}{Whether \code{weighted_neighbor_prediction} excludes zero-valued donors.}
\item{result}{A \code{DropoutKillerResult} with \code{uncertainty_available=TRUE}. Results from engines without predictive uncertainty are rejected.}
\item{n}{Number of uncertainty-aware completed matrices to draw.}
\item{seed}{Optional random seed for completed-matrix draws.}
}
\details{Factor-state learning uses a truncated SVD of standardized non-target genes. Every target gene carrying a recovery event is excluded from the factor-feature matrix, preventing direct target leakage. Ridge factor contribution is selected with exact analytic leave-one-out shrinkage. Under the default positive-conditional target, predictive variance is based on leave-one-out residual error and completed-matrix draws use a Gaussian approximation truncated below zero. Positive-conditional recovery is appropriate only for coordinates already classified as technical dropout; applying it to arbitrary natural zeros can bias expression upward.}
\value{\code{masked_factor_prediction} and \code{weighted_neighbor_prediction} return sparse predictions unless event tables are requested. \code{recover_dropout_expression} returns a selectively recovered matrix or a list with event diagnostics. \code{sample_dropout_expression} returns one completed matrix when \code{n=1} and a list otherwise.}
