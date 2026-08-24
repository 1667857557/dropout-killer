\name{local_alra_detect}
\alias{local_alra_detect}
\alias{local_alra_score}
\alias{select_dropout_mask}
\title{Membership-local low-rank dropout detection}
\description{Computes low-rank reconstructions independently within memberships. By default, biological-zero reconstruction variance is estimated from negative reconstructed values, shrunk toward a robust membership-level variance center, and used for one-sided zero-null tests with gene-wise Benjamini-Hochberg adjustment. The historical ALRA empirical-quantile gate remains available explicitly.}
\usage{
local_alra_detect(x, membership, rank = "auto", max_rank = 20L,
  rank_z = 6, quantile_prob = 0.001, min_cells = 8L,
  min_negative = 3L, seed = 1L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10)

local_alra_score(x, membership, rank = "auto", max_rank = 20L,
  rank_z = 6, quantile_prob = 0.001, min_cells = 8L,
  min_negative = 3L, seed = 1L,
  detection_method = c("eb_zero_null", "alra_quantile"),
  variance_prior_df = 10)

select_dropout_mask(score, threshold = 0.95, x = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix.}
\item{membership}{One membership label per cell.}
\item{rank}{Positive integer or \code{"auto"}.}
\item{max_rank}{Maximum rank considered in automatic selection.}
\item{rank_z}{Singular-spacing signal threshold in tail standard deviations for automatic rank selection.}
\item{quantile_prob}{Lower empirical quantile used only when \code{detection_method="alra_quantile"}. Retained for historical reproducibility.}
\item{min_cells}{Minimum membership size for local decomposition.}
\item{min_negative}{Minimum negative reconstructed values used to identify genes contributing to the robust membership-level zero-null variance center. Genes with fewer negatives can still be scored through variance shrinkage.}
\item{seed}{Random seed used by randomized SVD when required.}
\item{detection_method}{Detection engine. \code{"eb_zero_null"} is the default finite-sample detector; \code{"alra_quantile"} reproduces the historical empirical-tail gate.}
\item{variance_prior_df}{Non-negative equivalent prior sample size controlling shrinkage of each gene's negative-tail variance toward the membership-level robust center. Larger values increase shrinkage. Used only by \code{"eb_zero_null"}.}
\item{score}{Detection object or score matrix.}
\item{threshold}{Evidence cutoff in \code{(0,1]}. Under the default detector, event confidence is \code{1 - q_value}; therefore the default 0.95 selects approximately gene-wise BH q <= 0.05 under the working null.}
}
\details{For \code{detection_method="eb_zero_null"}, negative low-rank reconstruction values define a working symmetric biological-zero error distribution. Gene-specific second-moment variances are shrunk toward the median supported gene variance within the membership. Positive reconstructed observed zeros receive one-sided Gaussian p-values and are Benjamini-Hochberg adjusted within gene and membership, using the full number of observed zeros for that gene as the number of tests. The stored \code{confidence} equals \code{1-q_value}; it is not a calibrated Bayesian posterior probability of dropout. Gene-wise BH control is conditional on the working low-rank Gaussian zero-null and is not claimed to be exact global FDR control across all genes and memberships.}
\value{\code{local_alra_detect} returns a \code{DropoutKillerDetection} object containing event-level low-rank values, zero-null scale, z score, p value, q value, confidence, and membership diagnostics. \code{local_alra_score} returns a sparse confidence matrix. \code{select_dropout_mask} returns a sparse logical mask.}
