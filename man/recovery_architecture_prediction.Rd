\name{recovery_architecture_prediction}
\alias{recovery_architecture_prediction}
\title{Evaluate target-safe recovery architectures on a supplied dropout mask}
\usage{
recovery_architecture_prediction(
  x, counts, membership, embedding, mask, membership_fit = NULL,
  hard_stratum = NULL, methods = .dk_architecture_methods(),
  factor_rank = 5L, factor_features = 2000L, factor_ridge = 2,
  min_feature_observed = 20L, min_target_observed = 8L,
  factor_crossfit_folds = 5L, factor_crossfit_seed = 1L,
  bias_kappa = 10, architecture_geometry = NULL,
  return_details = TRUE
)
}
\arguments{
  \item{x}{Library-size normalized, log-transformed gene-by-cell expression.}
  \item{counts}{Raw UMI count matrix aligned to \code{x}.}
  \item{membership}{Membership vector or fitted membership assignment.}
  \item{embedding}{Cell-by-dimension biological embedding.}
  \item{mask}{Logical sparse matrix containing only supplied dropout targets.}
  \item{membership_fit}{Optional \code{DropoutKillerMembership} object.}
  \item{hard_stratum}{Optional hard biological stratum such as cell type.}
  \item{methods}{One or more architecture names.}
  \item{factor_rank, factor_features, factor_ridge}{P1 factor controls.}
  \item{min_feature_observed, min_target_observed}{Minimum support controls.}
  \item{factor_crossfit_folds, factor_crossfit_seed}{Target-fold cross-fitting controls.}
  \item{bias_kappa}{LOO bias-calibration shrinkage constant.}
  \item{architecture_geometry}{Optional reusable geometry state.}
  \item{return_details}{Return predictions, uncertainty, diagnostics, and state.}
}
\value{A list with event-level prediction and predictive-standard-deviation matrices.}
\description{
Builds the cross-fitted P1 state once and reuses it across local residual,
LOO blend, hierarchy partial-pooling, simplex, stabilized-state, fast P1, and
count generalized-Bayes recovery routes. Only supplied mask coordinates are
predicted; observed coordinates are not modified.
}
