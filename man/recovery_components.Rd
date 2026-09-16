\name{recover_dropout_expression}
\alias{recover_dropout_expression}
\alias{sample_dropout_expression}
\title{Selective dropout recovery and uncertainty sampling}
\description{Recovers only supplied zero-mask events without overwriting observed coordinates. P1 stabilized state is the production default; comparison recovery engines remain selectable through this expert interface for controlled ablation.}
\usage{
recover_dropout_expression(x, mask, membership, embedding = NULL,
  neighbor_k = 30L, neighbor_sigma = NULL,
  min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
  cap_quantile = NULL, return_details = FALSE,
  recovery_method = c("p1_stabilized_state", "masked_factor",
  "tree_local_factor", "neighbor"), factor_rank = 5L,
  factor_features = 2000L, factor_ridge = 2,
  min_feature_observed = 20L, min_target_observed = 8L,
  factor_target = c("positive", "all_observed"),
  membership_fit = NULL, hard_stratum = NULL,
  tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
  candidate_k = 100L, min_effective_donors = 5,
  local_info_kappa = 5, factor_crossfit_folds = 5L,
  factor_crossfit_seed = 1L, support_adaptive_rank = TRUE,
  bias_kappa = 10, predictor_smoothing = 0.25)

sample_dropout_expression(result, n = 1L, seed = NULL)
}
\arguments{
\item{x}{Non-negative gene-by-cell expression matrix on the desired recovery scale.}
\item{mask}{Sparse or dense logical matrix containing only zero-valued recovery events.}
\item{membership}{One membership label per cell, or a \code{DropoutKillerMembership} object.}
\item{embedding}{Cell-by-dimension biological representation. Required for P1 stabilized state, neighbor, and tree-local recovery.}
\item{neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only}{Controls for the optional neighbor comparison engine.}
\item{cap_quantile}{Optional upper cap used by comparison recovery engines.}
\item{return_details}{Return event diagnostics and local geometry with the recovered matrix.}
\item{recovery_method}{Recovery engine. \code{"p1_stabilized_state"} is the production default. Other values are retained for controlled recovery ablation.}
\item{factor_rank, factor_features, factor_ridge}{Target-safe coexpression-state and ridge settings.}
\item{min_feature_observed}{Minimum observations for a predictor feature.}
\item{min_target_observed}{Minimum positive target donors required before a factor state is attempted.}
\item{factor_target}{Recovery target estimand. P1 stabilized state and tree-local recovery support the positive-conditional target.}
\item{membership_fit}{Optional retained hierarchy associated with \code{membership}.}
\item{hard_stratum}{Optional absolute biological borrowing boundary for direct recovery.}
\item{tree_weight, tree_tau, local_k, candidate_k, min_effective_donors, local_info_kappa}{Expert controls for hierarchy-aware comparison recovery.}
\item{factor_crossfit_folds}{Number of deterministic target-gene folds.}
\item{factor_crossfit_seed}{Seed controlling target-gene fold assignment.}
\item{support_adaptive_rank}{Cap target regression rank using available positive-donor support.}
\item{bias_kappa}{Leave-one-out bias-shrinkage constant.}
\item{predictor_smoothing}{Weight in \eqn{[0,1]} for one row-stochastic smoothing step applied to the target-safe predictor state.}
\item{result}{A \code{DropoutKillerResult} with predictive uncertainty.}
\item{n}{Number of uncertainty-aware completed matrices.}
\item{seed}{Optional random seed.}
}
\details{P1 stabilized-state recovery excludes each target-gene fold before predictor construction, smooths only the target-safe predictor state over biological geometry, learns a low-dimensional cell state, and fits positive-donor ridge predictions with analytic leave-one-out shrinkage and bias calibration. Observed expression is never smoothed or overwritten. Positive completed-matrix draws use Gamma moment matching when predictive variance is available.}
\value{\code{recover_dropout_expression} returns a selectively recovered matrix or, with \code{return_details=TRUE}, a list containing expression, event diagnostics, uncertainty status, and local geometry. \code{sample_dropout_expression} returns one or more completed matrices.}
