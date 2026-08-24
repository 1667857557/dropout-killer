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
       n_donors = cp$n_donors, bandwidth = cp$bandwidth)
}

.dk_recover_events <- function(x, events, membership, embedding = NULL,
                               recovery_method = c("masked_factor", "neighbor"),
                               factor_rank = 5L, factor_features = 2000L,
                               factor_ridge = 1, min_feature_observed = 20L,
                               min_target_observed = 20L,
                               neighbor_k = 30L, neighbor_sigma = NULL,
                               min_positive_neighbors = 1L,
                               neighbor_positive_only = TRUE,
                               cap_quantile = NULL,
                               factor_target = c("positive", "all_observed")) {
  recovery_method <- match.arg(recovery_method)
  factor_target <- match.arg(factor_target)
  if (recovery_method == "masked_factor") {
    mf <- .dk_masked_factor_predict_events(
      x, membership, events, factor_rank = factor_rank,
      factor_features = factor_features, factor_ridge = factor_ridge,
      min_feature_observed = min_feature_observed,
      min_target_observed = min_target_observed,
      cap_quantile = cap_quantile, target_mode = factor_target
    )
    n <- nrow(events)
    return(c(mf, list(cell_prediction = rep(NA_real_, n),
                      cell_available = rep(FALSE, n),
                      n_donors = integer(n), bandwidth = rep(NA_real_, n))))
  }
  if (is.null(embedding)) stop("embedding is required for recovery_method='neighbor'", call. = FALSE)
  .dk_neighbor_recover_events(
    x, events, membership, embedding, neighbor_k, neighbor_sigma,
    min_positive_neighbors, neighbor_positive_only, cap_quantile
  )
}

#' Selectively recover masked zero events
#'
#' The default masked-factor engine treats supplied dropout events as missing,
#' learns membership-local coexpression factors from non-target genes, and
#' predicts target magnitude from reliable positive donors. Thus, once a zero
#' has already been classified as a technical dropout, recovery estimates
#' `E[X_g | X_g > 0, cell state]` rather than the zero-inflated unconditional
#' membership mean. Set `factor_target = "all_observed"` to reproduce the
#' previous target. Observed values are never overwritten.
#'
#' @export
recover_dropout_expression <- function(x, mask, membership, embedding = NULL,
                                       neighbor_k = 30L, neighbor_sigma = NULL,
                                       min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                                       cap_quantile = NULL, return_details = FALSE,
                                       recovery_method = c("masked_factor", "neighbor"),
                                       factor_rank = 5L, factor_features = 2000L,
                                       factor_ridge = 1, min_feature_observed = 20L,
                                       min_target_observed = 20L,
                                       factor_target = c("positive", "all_observed")) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  recovery_method <- match.arg(recovery_method)
  factor_target <- match.arg(factor_target)
  z <- if (recovery_method == "neighbor") .dk_align_embedding(embedding, nm$cells) else NULL
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
    cap_quantile = cap_quantile, factor_target = factor_target
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
  events$recovered <- fit$prediction
  events$changed <- ok
  uncertainty_available <- identical(recovery_method, "masked_factor") &&
    (!any(ok) || all(is.finite(fit$prediction_sd[ok]) & fit$prediction_sd[ok] >= 0))
  if (return_details) {
    list(expression = out, events = events,
         uncertainty_available = uncertainty_available,
         factor_target = factor_target)
  } else out
}

#' Run selective dropout detection and recovery
#'
#' The default detector keeps membership-local low-rank structure but replaces
#' the finite-sample-unstable empirical extreme-quantile gate with an
#' empirical-Bayes symmetric zero-null test. Gene-specific negative-tail
#' variance is shrunk toward a robust membership center and positive
#' reconstructed zeros are tested with gene-wise BH correction.
#'
#' Recovery then conditions on the selected event being a technical dropout and
#' predicts positive expression magnitude from membership-local coexpression
#' factors. Observed values, including unmasked biological zeros, are never
#' overwritten.
#'
#' @export
dropout_killer <- function(x, embedding, membership = NULL, group = NULL, split_by = NULL,
                           gamma = 150, k_knn = 5L, approximate = "auto", approx_n = 20000L,
                           rank = "auto", max_rank = 20L, rank_z = 6, quantile_prob = 0.001,
                           threshold = 0.95, min_cells = 8L, min_negative = 3L,
                           neighbor_k = 30L, neighbor_sigma = NULL,
                           min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                           cap_quantile = NULL, seed = 12345L, return_score = FALSE,
                           recovery_method = c("masked_factor", "neighbor"),
                           factor_rank = 5L, factor_features = 2000L,
                           factor_ridge = 1, min_feature_observed = 20L,
                           min_target_observed = 20L,
                           detection_method = c("eb_zero_null", "alra_quantile"),
                           variance_prior_df = 10,
                           factor_target = c("positive", "all_observed")) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  z <- .dk_align_embedding(embedding, nm$cells)
  recovery_method <- match.arg(recovery_method)
  detection_method <- match.arg(detection_method)
  factor_target <- match.arg(factor_target)
  group <- .dk_align_vector(group, nm$cells, "group")
  split_by <- .dk_align_vector(split_by, nm$cells, "split_by")
  membership_fit <- NULL
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
    x, ev, membership, if (recovery_method == "neighbor") z else NULL,
    recovery_method = recovery_method,
    factor_rank = factor_rank, factor_features = factor_features,
    factor_ridge = factor_ridge, min_feature_observed = min_feature_observed,
    min_target_observed = min_target_observed,
    neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
    min_positive_neighbors = min_positive_neighbors,
    neighbor_positive_only = neighbor_positive_only,
    cap_quantile = cap_quantile, factor_target = factor_target
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
    events$bandwidth <- numeric(); events$recovered <- numeric()
    events$changed <- logical()
  }

  uncertainty_available <- identical(recovery_method, "masked_factor") &&
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
    mask = mask, events = events, predictive_variance = predictive_variance,
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
  if (!is.null(x$settings$factor_target) && x$settings$recovery_method == "masked_factor")
    cat(" recovery target:", x$settings$factor_target, "\n")
  cat(" high-confidence dropout events:", nrow(x$events), "\n")
  cat(" recovered events:", if (nrow(x$events)) sum(x$events$changed) else 0L, "\n")
  invisible(x)
}
