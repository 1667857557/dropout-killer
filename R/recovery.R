.dk_neighbor_recover_events <- function(x, events, membership, embedding,
                                        neighbor_k, neighbor_sigma, min_positive_neighbors,
                                        neighbor_positive_only, cap_quantile = NULL) {
  if (!is.logical(neighbor_positive_only) || length(neighbor_positive_only) != 1L || is.na(neighbor_positive_only)) stop("neighbor_positive_only must be TRUE or FALSE", call. = FALSE)
  cp <- .dk_cell_predict_events(x, embedding, membership, events, neighbor_k, neighbor_sigma,
                                min_positive_neighbors, positive_only = neighbor_positive_only)
  available <- is.finite(cp$prediction); pred <- cp$prediction
  pred[!available | pred < 0] <- 0
  if (!is.null(cap_quantile)) {
    if (!is.numeric(cap_quantile) || length(cap_quantile) != 1L || !is.finite(cap_quantile) || cap_quantile <= 0.5 || cap_quantile > 1) stop("cap_quantile must be in (0.5,1]", call. = FALSE)
    for (m in unique(events$membership)) {
      q <- which(events$membership == m); cells <- which(membership == m)
      for (g in unique(events$i[q])) {
        z <- q[events$i[q] == g]; obs <- as.numeric(x[g, cells]); obs <- obs[obs > 0]
        if (length(obs)) pred[z] <- pmin(pred[z], as.numeric(stats::quantile(obs, cap_quantile, names = FALSE, type = 8)))
      }
    }
  }
  n <- nrow(events)
  list(prediction = pred, factor_prediction = rep(NA_real_, n), prediction_sd = rep(NA_real_, n),
       predictability = rep(NA_real_, n), shrinkage = rep(NA_real_, n),
       n_observed_gene = integer(n), factor_rank = integer(n), factor_features = integer(n),
       factor_iterations = integer(n), factor_converged = logical(n),
       recovery_method = rep("neighbor", n), target_mode = rep("neighbor", n),
       cell_prediction = cp$prediction, cell_available = available,
       n_donors = cp$n_donors, bandwidth = cp$bandwidth,
       local_positive_mean = rep(NA_real_, n), local_positive_variance = rep(NA_real_, n),
       local_positive_prevalence = rep(NA_real_, n), effective_donors = rep(NA_real_, n),
       tree_distance_weighted_mean = rep(NA_real_, n),
       embedding_distance_weighted_mean = rep(NA_real_, n), geometry = NULL)
}

.dk_recover_events <- function(x, events, membership, embedding = NULL,
                               recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
                               factor_rank = 5L, factor_features = 2000L,
                               factor_ridge = 1, min_feature_observed = 20L,
                               min_target_observed = 20L,
                               neighbor_k = 30L, neighbor_sigma = NULL,
                               min_positive_neighbors = 1L,
                               neighbor_positive_only = TRUE,
                               cap_quantile = NULL,
                               factor_target = c("positive", "all_observed"),
                               membership_fit = NULL, hard_stratum = NULL,
                               tree_weight = 0.5, tree_tau = NULL,
                               local_k = 30L, candidate_k = 100L,
                               min_effective_donors = 5,
                               local_info_kappa = 5) {
  recovery_method <- match.arg(recovery_method)
  factor_target <- match.arg(factor_target)
  n <- nrow(events)
  if (recovery_method == "masked_factor") {
    mf <- .dk_masked_factor_predict_events(
      x, membership, events, factor_rank = factor_rank,
      factor_features = factor_features, factor_ridge = factor_ridge,
      min_feature_observed = min_feature_observed,
      min_target_observed = min_target_observed,
      cap_quantile = cap_quantile, target_mode = factor_target
    )
    return(c(mf, list(cell_prediction = rep(NA_real_, n),
                      cell_available = rep(FALSE, n), n_donors = integer(n),
                      bandwidth = rep(NA_real_, n), local_positive_mean = rep(NA_real_, n),
                      local_positive_variance = rep(NA_real_, n),
                      local_positive_prevalence = rep(NA_real_, n),
                      effective_donors = rep(NA_real_, n),
                      tree_distance_weighted_mean = rep(NA_real_, n),
                      embedding_distance_weighted_mean = rep(NA_real_, n), geometry = NULL)))
  }
  if (recovery_method == "tree_local_factor") {
    if (factor_target != "positive") stop("tree_local_factor implements the positive-conditional recovery target only", call. = FALSE)
    if (is.null(embedding)) stop("embedding is required for recovery_method='tree_local_factor'", call. = FALSE)
    tl <- .dk_tree_local_predict_events(
      x, embedding, membership, events, membership_fit = membership_fit,
      hard_stratum = hard_stratum, factor_rank = factor_rank,
      factor_features = factor_features, factor_ridge = factor_ridge,
      min_feature_observed = min_feature_observed,
      min_target_observed = min_target_observed,
      tree_weight = tree_weight, tree_tau = tree_tau,
      local_k = local_k, candidate_k = candidate_k,
      min_effective_donors = min_effective_donors,
      local_info_kappa = local_info_kappa
    )
    tl$cell_prediction <- tl$local_positive_mean
    tl$cell_available <- is.finite(tl$local_positive_mean)
    tl$n_donors <- tl$n_observed_gene
    tl$bandwidth <- if (!is.null(tl$geometry)) tl$geometry$bandwidth[events$j] else rep(NA_real_, n)
    return(tl)
  }
  if (is.null(embedding)) stop("embedding is required for recovery_method='neighbor'", call. = FALSE)
  .dk_neighbor_recover_events(
    x, events, membership, embedding, neighbor_k, neighbor_sigma,
    min_positive_neighbors, neighbor_positive_only, cap_quantile
  )
}

#' Selectively recover masked zero events
#'
#' `tree_local_factor` preserves the SuperCell biological geometry instead of
#' treating every cell in a final membership as exchangeable. Hard biological
#' strata remain absolute borrowing boundaries; within a stratum, the retained
#' walktrap hierarchy and the original embedding assign larger weights to
#' biologically closer donors. The positive target baseline is therefore
#' query-specific. Coexpression factors predict only residual expression beyond
#' that local baseline. `masked_factor` and `neighbor` remain available as
#' comparison engines. Observed values are never overwritten.
#'
#' @export
recover_dropout_expression <- function(x, mask, membership, embedding = NULL,
                                       neighbor_k = 30L, neighbor_sigma = NULL,
                                       min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                                       cap_quantile = NULL, return_details = FALSE,
                                       recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
                                       factor_rank = 5L, factor_features = 2000L,
                                       factor_ridge = 1, min_feature_observed = 20L,
                                       min_target_observed = 20L,
                                       factor_target = c("positive", "all_observed"),
                                       membership_fit = NULL, hard_stratum = NULL,
                                       tree_weight = 0.5, tree_tau = NULL,
                                       local_k = 30L, candidate_k = 100L,
                                       min_effective_donors = 5,
                                       local_info_kappa = 5) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  if (inherits(membership, "DropoutKillerMembership")) {
    membership_fit <- membership; membership <- membership_fit$membership
  }
  membership <- .dk_align_membership(membership, nm$cells)
  recovery_method <- match.arg(recovery_method)
  factor_target <- match.arg(factor_target)
  z <- if (recovery_method %in% c("neighbor", "tree_local_factor")) .dk_align_embedding(embedding, nm$cells) else NULL
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x)))) stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) {
    if (any(events$i > nrow(x) | events$j > ncol(x))) stop("mask dimensions exceed x", call. = FALSE)
    is_zero <- as.vector(x[cbind(events$i, events$j)] == 0)
    if (!all(is_zero)) stop("mask contains observed non-zero entries; selective recovery refuses to overwrite them", call. = FALSE)
    events$membership <- membership[events$j]
  } else events$membership <- integer()
  fit <- .dk_recover_events(
    x, events, membership, z, recovery_method = recovery_method,
    factor_rank = factor_rank, factor_features = factor_features,
    factor_ridge = factor_ridge, min_feature_observed = min_feature_observed,
    min_target_observed = min_target_observed,
    neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
    min_positive_neighbors = min_positive_neighbors,
    neighbor_positive_only = neighbor_positive_only,
    cap_quantile = cap_quantile, factor_target = factor_target,
    membership_fit = membership_fit, hard_stratum = hard_stratum,
    tree_weight = tree_weight, tree_tau = tree_tau,
    local_k = local_k, candidate_k = candidate_k,
    min_effective_donors = min_effective_donors,
    local_info_kappa = local_info_kappa
  )
  ok <- fit$prediction > 0 & is.finite(fit$prediction)
  if (inherits(x, "Matrix")) {
    delta <- .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok],
                                nrow(x), ncol(x), dimnames(x))
    out <- x + delta
  } else {
    out <- x
    if (any(ok)) out[cbind(events$i[ok], events$j[ok])] <- fit$prediction[ok]
  }
  events$factor_prediction <- fit$factor_prediction
  events$prediction_sd <- fit$prediction_sd
  events$predictability <- fit$predictability
  events$shrinkage <- fit$shrinkage
  events$n_observed_gene <- fit$n_observed_gene
  events$factor_rank <- fit$factor_rank
  events$factor_features <- fit$factor_features
  events$factor_iterations <- fit$factor_iterations
  events$factor_converged <- fit$factor_converged
  events$recovery_method <- fit$recovery_method
  events$target_mode <- fit$target_mode
  events$cell_prediction <- fit$cell_prediction
  events$cell_available <- fit$cell_available
  events$n_donors <- fit$n_donors
  events$bandwidth <- fit$bandwidth
  events$local_positive_mean <- fit$local_positive_mean
  events$local_positive_variance <- fit$local_positive_variance
  events$local_positive_prevalence <- fit$local_positive_prevalence
  events$effective_donors <- fit$effective_donors
  events$tree_distance_weighted_mean <- fit$tree_distance_weighted_mean
  events$embedding_distance_weighted_mean <- fit$embedding_distance_weighted_mean
  events$recovered <- fit$prediction
  events$changed <- ok
  uncertainty_available <- recovery_method %in% c("masked_factor", "tree_local_factor") &&
    (!any(ok) || all(is.finite(fit$prediction_sd[ok]) & fit$prediction_sd[ok] >= 0))
  if (return_details) {
    list(expression = out, events = events, uncertainty_available = uncertainty_available,
         factor_target = factor_target, local_geometry = fit$geometry)
  } else out
}

#' Run selective dropout detection and recovery
#'
#' Detection remains separate from recovery. The tree-local engine retains the
#' full SuperCell walktrap hierarchy generated inside each hard biological
#' stratum. A final gamma cut defines a high-weight core, not an absolute wall:
#' nearby sibling memberships can contribute with rapidly decaying hierarchy and
#' embedding weights, while different hard strata never borrow from one another.
#'
#' @export
dropout_killer <- function(x, embedding, membership = NULL, group = NULL, split_by = NULL,
                           gamma = 150, k_knn = 5L, approximate = "auto", approx_n = 20000L,
                           rank = "auto", max_rank = 20L, rank_z = 6, quantile_prob = 0.001,
                           threshold = 0.95, min_cells = 8L, min_negative = 3L,
                           neighbor_k = 30L, neighbor_sigma = NULL,
                           min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                           cap_quantile = NULL, seed = 12345L, return_score = FALSE,
                           recovery_method = c("masked_factor", "tree_local_factor", "neighbor"),
                           factor_rank = 5L, factor_features = 2000L,
                           factor_ridge = 1, min_feature_observed = 20L,
                           min_target_observed = 20L,
                           detection_method = c("eb_zero_null", "alra_quantile"),
                           variance_prior_df = 10,
                           factor_target = c("positive", "all_observed"),
                           tree_weight = 0.5, tree_tau = NULL,
                           local_k = 30L, candidate_k = 100L,
                           min_effective_donors = 5,
                           local_info_kappa = 5) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  z <- .dk_align_embedding(embedding, nm$cells)
  recovery_method <- match.arg(recovery_method)
  detection_method <- match.arg(detection_method)
  factor_target <- match.arg(factor_target)
  group <- .dk_align_vector(group, nm$cells, "group")
  split_by <- .dk_align_vector(split_by, nm$cells, "split_by")
  hard_recovery_stratum <- if (!is.null(group) || !is.null(split_by)) .dk_stratum(group, split_by, ncol(x)) else NULL
  membership_fit <- NULL
  if (inherits(membership, "DropoutKillerMembership")) {
    membership_fit <- membership; membership <- membership_fit$membership
  }
  if (is.null(membership)) {
    membership_fit <- build_supercell_membership(
      z, group = group, split_by = split_by, gamma = gamma,
      k_knn = k_knn, approximate = approximate,
      approx_n = approx_n, seed = seed
    )
    membership <- .dk_align_membership(membership_fit$membership, nm$cells)
    membership_fit$membership <- membership
    membership_fit$membership_table <- as.data.frame(table(membership), stringsAsFactors = FALSE)
    names(membership_fit$membership_table) <- c("membership", "n_cells")
    membership_fit$membership_table$membership <- as.integer(as.character(membership_fit$membership_table$membership))
  } else membership <- .dk_align_membership(membership, nm$cells)

  det <- local_alra_detect(
    x, membership, rank = rank, max_rank = max_rank, rank_z = rank_z,
    quantile_prob = quantile_prob, min_cells = min_cells,
    min_negative = min_negative, seed = seed,
    detection_method = detection_method,
    variance_prior_df = variance_prior_df
  )
  mask <- select_dropout_mask(det, threshold = threshold)
  ev <- .dk_mask_events(mask)
  if (nrow(ev)) ev$membership <- membership[ev$j] else ev$membership <- integer()

  rec <- .dk_recover_events(
    x, ev, membership, if (recovery_method %in% c("neighbor", "tree_local_factor")) z else NULL,
    recovery_method = recovery_method,
    factor_rank = factor_rank, factor_features = factor_features,
    factor_ridge = factor_ridge, min_feature_observed = min_feature_observed,
    min_target_observed = min_target_observed,
    neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
    min_positive_neighbors = min_positive_neighbors,
    neighbor_positive_only = neighbor_positive_only,
    cap_quantile = cap_quantile, factor_target = factor_target,
    membership_fit = membership_fit, hard_stratum = hard_recovery_stratum,
    tree_weight = tree_weight, tree_tau = tree_tau,
    local_k = local_k, candidate_k = candidate_k,
    min_effective_donors = min_effective_donors,
    local_info_kappa = local_info_kappa
  )
  ok <- rec$prediction > 0 & is.finite(rec$prediction)
  if (inherits(x, "Matrix")) {
    delta <- .dk_sparse_numeric(ev$i[ok], ev$j[ok], rec$prediction[ok],
                                nrow(x), ncol(x), dimnames(x))
    expression <- x + delta
  } else {
    expression <- x
    if (any(ok)) expression[cbind(ev$i[ok], ev$j[ok])] <- rec$prediction[ok]
  }

  if (nrow(ev)) {
    key <- paste(ev$i, ev$j, sep = ":")
    did <- match(key, paste(det$events$i, det$events$j, sep = ":"))
    events <- det$events[did, , drop = FALSE]
    events$factor_prediction <- rec$factor_prediction
    events$prediction_sd <- rec$prediction_sd
    events$predictability <- rec$predictability
    events$shrinkage <- rec$shrinkage
    events$n_observed_gene <- rec$n_observed_gene
    events$factor_rank <- rec$factor_rank
    events$factor_features <- rec$factor_features
    events$factor_iterations <- rec$factor_iterations
    events$factor_converged <- rec$factor_converged
    events$recovery_method <- rec$recovery_method
    events$target_mode <- rec$target_mode
    events$cell_prediction <- rec$cell_prediction
    events$cell_available <- rec$cell_available
    events$n_donors <- rec$n_donors
    events$bandwidth <- rec$bandwidth
    events$local_positive_mean <- rec$local_positive_mean
    events$local_positive_variance <- rec$local_positive_variance
    events$local_positive_prevalence <- rec$local_positive_prevalence
    events$effective_donors <- rec$effective_donors
    events$tree_distance_weighted_mean <- rec$tree_distance_weighted_mean
    events$embedding_distance_weighted_mean <- rec$embedding_distance_weighted_mean
    events$recovered <- rec$prediction
    events$changed <- ok
  } else {
    events <- det$events[FALSE, , drop = FALSE]
    events$factor_prediction <- numeric(); events$prediction_sd <- numeric()
    events$predictability <- numeric(); events$shrinkage <- numeric()
    events$n_observed_gene <- integer(); events$factor_rank <- integer()
    events$factor_features <- integer(); events$factor_iterations <- integer()
    events$factor_converged <- logical(); events$recovery_method <- character()
    events$target_mode <- character(); events$cell_prediction <- numeric()
    events$cell_available <- logical(); events$n_donors <- integer()
    events$bandwidth <- numeric(); events$local_positive_mean <- numeric()
    events$local_positive_variance <- numeric(); events$local_positive_prevalence <- numeric()
    events$effective_donors <- numeric(); events$tree_distance_weighted_mean <- numeric()
    events$embedding_distance_weighted_mean <- numeric(); events$recovered <- numeric()
    events$changed <- logical()
  }

  uncertainty_available <- recovery_method %in% c("masked_factor", "tree_local_factor") &&
    (!any(ok) || all(is.finite(rec$prediction_sd[ok]) & rec$prediction_sd[ok] >= 0))
  predictive_variance <- NULL
  if (uncertainty_available) {
    pv_ok <- ok & is.finite(rec$prediction_sd) & rec$prediction_sd >= 0
    predictive_variance <- .dk_sparse_numeric(
      ev$i[pv_ok], ev$j[pv_ok], rec$prediction_sd[pv_ok]^2,
      nrow(x), ncol(x), dimnames(x)
    )
  }

  out <- list(
    expression = expression, membership = membership, membership_fit = membership_fit,
    local_geometry = rec$geometry, mask = mask, events = events,
    predictive_variance = predictive_variance,
    uncertainty_available = uncertainty_available, detection = det,
    settings = list(
      gamma = gamma, k_knn = k_knn, approximate = approximate, approx_n = approx_n,
      rank = rank, max_rank = max_rank, rank_z = rank_z,
      quantile_prob = quantile_prob, threshold = threshold,
      min_cells = min_cells, min_negative = min_negative,
      detection_method = detection_method, variance_prior_df = variance_prior_df,
      recovery_method = recovery_method, factor_target = factor_target,
      factor_rank = factor_rank, factor_features = factor_features,
      factor_ridge = factor_ridge,
      min_feature_observed = min_feature_observed,
      min_target_observed = min_target_observed,
      tree_weight = tree_weight, tree_tau = tree_tau,
      local_k = local_k, candidate_k = candidate_k,
      min_effective_donors = min_effective_donors,
      local_info_kappa = local_info_kappa,
      neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
      min_positive_neighbors = min_positive_neighbors,
      neighbor_positive_only = neighbor_positive_only,
      cap_quantile = cap_quantile, seed = seed
    )
  )
  if (return_score) {
    sev <- det$events
    out$score <- .dk_sparse_numeric(
      sev$i, sev$j, sev$confidence, det$dimensions[1L],
      det$dimensions[2L], det$dimnames
    )
    attr(out$score, "zero_only") <- TRUE
    attr(out$score, "detection") <- det
  }
  class(out) <- "DropoutKillerResult"
  out
}

#' @export
DropoutKiller <- dropout_killer
#' @export
run_dropout_killer <- dropout_killer

#' @export
print.DropoutKillerResult <- function(x, ...) {
  cat("DropoutKiller result\n")
  cat(" dimensions:", paste(dim(x$expression), collapse = " x "), "\n")
  cat(" memberships:", length(unique(x$membership)), "\n")
  cat(" detector:", x$settings$detection_method, "\n")
  cat(" recovery engine:", x$settings$recovery_method, "\n")
  if (!is.null(x$settings$factor_target) && x$settings$recovery_method %in% c("masked_factor", "tree_local_factor"))
    cat(" recovery target:", x$settings$factor_target, "\n")
  cat(" high-confidence dropout events:", nrow(x$events), "\n")
  cat(" recovered events:", if (nrow(x$events)) sum(x$events$changed) else 0L, "\n")
  invisible(x)
}
