#' Run selective dropout detection and recovery
#'
#' The default zero detector is a source-faithful R port of the scGACL
#' Gamma-Normal mixture detector. It operates on the supplied raw count matrix,
#' transforms counts as `log(1.01 + count)`, fits one Gamma-Normal mixture for
#' every gene within each supplied `group` (the released scGACL default uses
#' known cell-type labels), and selects observed zeros whose Gamma-component
#' posterior is at least `scgacl_dropout_threshold = 0.5`.
#'
#' Detection and recovery intentionally use separate input scales. scGACL sees
#' the original raw counts. Recovery continues to use the existing ALRA
#' library-size-to-10,000 plus `log1p` working matrix when `normalize = TRUE`.
#'
#' `detection_method = "alra_global"` is the optional ALRA comparator. It follows
#' the released KlugerLab/ALRA implementation as one all-cell matrix: original
#' ALRA normalization, singular-value-spacing rank selection, randomized SVD,
#' and the gene-wise absolute 0.1% low-rank quantile gate. Cell classes and
#' SuperCell memberships do not fragment this ALRA comparator.
#'
#' The production recovery engine remains `p1_stabilized_state`.
#'
#' @export
dropout_killer <- function(x, embedding, membership = NULL, group = NULL, split_by = NULL,
                           gamma = 150, k_knn = 5L, approximate = "auto", approx_n = 20000L,
                           rank = "auto", max_rank = 20L, rank_z = 6, quantile_prob = 0.001,
                           threshold = 0.95, min_cells = 8L, min_negative = 3L,
                           neighbor_k = 30L, neighbor_sigma = NULL,
                           min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                           cap_quantile = NULL, seed = 12345L, return_score = FALSE,
                           recovery_method = c("p1_stabilized_state", "masked_factor",
                                               "tree_local_factor", "neighbor"),
                           factor_rank = 5L, factor_features = 2000L,
                           factor_ridge = 2, min_feature_observed = 20L,
                           min_target_observed = 8L,
                           detection_method = c("scgacl_gamma_normal", "alra_global",
                                                "alra_global_by_group", "eb_zero_null",
                                                "alra_quantile"),
                           scgacl_dropout_threshold = 0.5,
                           variance_prior_df = 10,
                           factor_target = c("positive", "all_observed"),
                           tree_weight = 0.5, tree_tau = NULL,
                           local_k = 30L, candidate_k = 100L,
                           min_effective_donors = 5,
                           local_info_kappa = 5,
                           normalize = TRUE, normalization_scale_factor = 1e4,
                           alra_K = 100L, alra_noise_start = 80L,
                           alra_choose_q = 2L, alra_svd_q = 10L,
                           factor_crossfit_folds = 5L,
                           factor_crossfit_seed = 1L,
                           support_adaptive_rank = TRUE,
                           bias_kappa = 10,
                           predictor_smoothing = 0.25) {
  x_raw <- .dk_validate_expression(x)
  nm <- .dk_names(x_raw)

  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("normalize must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(normalization_scale_factor) || length(normalization_scale_factor) != 1L ||
      !is.finite(normalization_scale_factor) || normalization_scale_factor <= 0) {
    stop("normalization_scale_factor must be a finite value > 0", call. = FALSE)
  }
  if (!is.numeric(scgacl_dropout_threshold) || length(scgacl_dropout_threshold) != 1L ||
      !is.finite(scgacl_dropout_threshold) || scgacl_dropout_threshold < 0 ||
      scgacl_dropout_threshold > 1) {
    stop("scgacl_dropout_threshold must be a finite value in [0,1]", call. = FALSE)
  }

  recovery_method <- match.arg(recovery_method)
  detection_method <- match.arg(detection_method)
  factor_target <- match.arg(factor_target)
  group <- .dk_align_vector(group, nm$cells, "group")
  split_by <- .dk_align_vector(split_by, nm$cells, "split_by")

  if (identical(detection_method, "alra_global_by_group")) {
    warning(
      "detection_method='alra_global_by_group' is deprecated. ALRA now follows ",
      "the original all-cell global method; use detection_method='alra_global'.",
      call. = FALSE
    )
    detection_method <- "alra_global"
  }

  if (identical(detection_method, "scgacl_gamma_normal")) {
    vals <- if (inherits(x_raw, "sparseMatrix") && "x" %in% methods::slotNames(x_raw)) {
      methods::slot(x_raw, "x")
    } else if (inherits(x_raw, "sparseMatrix")) {
      rep.int(1, nrow(Matrix::summary(x_raw)))
    } else {
      as.vector(x_raw)
    }
    if (length(vals) && any(abs(vals - round(vals)) > 1e-8)) {
      stop(
        "detection_method='scgacl_gamma_normal' requires raw count values. ",
        "The released scGACL detector applies log(1.01 + raw_count) internally.",
        call. = FALSE
      )
    }
  }

  # Recovery scale is independent from the scGACL detector scale.
  x_work <- if (normalize) {
    .dk_alra_library_log(x_raw, scale_factor = normalization_scale_factor)
  } else {
    x_raw
  }

  z <- .dk_align_embedding(embedding, nm$cells)
  hard_recovery_stratum <- if (!is.null(group) || !is.null(split_by)) {
    .dk_stratum(group, split_by, ncol(x_work))
  } else NULL

  membership_fit <- NULL
  if (inherits(membership, "DropoutKillerMembership")) {
    membership_fit <- membership
    membership <- membership_fit$membership
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
    membership_fit$membership_table$membership <- as.integer(
      as.character(membership_fit$membership_table$membership)
    )
  } else {
    membership <- .dk_align_membership(membership, nm$cells)
  }

  if (identical(detection_method, "scgacl_gamma_normal")) {
    det <- .dk_scgacl_detect(
      x_raw, group = group,
      dropout_threshold = scgacl_dropout_threshold,
      point = .dk_scgacl_point()
    )
    if (nrow(det$events)) det$events$membership <- membership[det$events$j]
    mask <- .dk_sparse_logical(
      det$events$i, det$events$j,
      nrow(x_work), ncol(x_work), dimnames(x_work)
    )
  } else if (identical(detection_method, "alra_global")) {
    # Original ALRA detection operates on its original library/log scale.
    # When normalize=TRUE x_work is exactly that scale. If normalize=FALSE, the
    # caller explicitly assumes responsibility for supplying an ALRA-compatible
    # working matrix, preserving the package's historical pre-normalized path.
    det <- .dk_original_alra_detect(
      x_work, rank = rank, quantile_prob = quantile_prob,
      seed = seed, K = alra_K, rank_z = rank_z,
      noise_start = alra_noise_start, choose_q = alra_choose_q,
      svd_q = alra_svd_q
    )
    if (nrow(det$events)) det$events$membership <- membership[det$events$j]
    mask <- .dk_sparse_logical(
      det$events$i, det$events$j,
      nrow(x_work), ncol(x_work), dimnames(x_work)
    )
  } else {
    det <- local_alra_detect(
      x_work, membership, rank = rank, max_rank = max_rank, rank_z = rank_z,
      quantile_prob = quantile_prob, min_cells = min_cells,
      min_negative = min_negative, seed = seed,
      detection_method = detection_method,
      variance_prior_df = variance_prior_df
    )
    mask <- select_dropout_mask(det, threshold = threshold)
  }

  ev <- .dk_mask_events(mask)
  if (nrow(ev)) ev$membership <- membership[ev$j] else ev$membership <- integer()

  rec <- .dk_recover_events(
    x_work, ev, membership, if (recovery_method %in% c(
      "p1_stabilized_state", "neighbor", "tree_local_factor"
    )) z else NULL,
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
    local_info_kappa = local_info_kappa,
    factor_crossfit_folds = factor_crossfit_folds,
    factor_crossfit_seed = factor_crossfit_seed,
    support_adaptive_rank = support_adaptive_rank,
    bias_kappa = bias_kappa,
    predictor_smoothing = predictor_smoothing
  )

  ok <- rec$prediction > 0 & is.finite(rec$prediction)
  if (inherits(x_work, "Matrix")) {
    delta <- .dk_sparse_numeric(
      ev$i[ok], ev$j[ok], rec$prediction[ok],
      nrow(x_work), ncol(x_work), dimnames(x_work)
    )
    expression <- x_work + delta
  } else {
    expression <- x_work
    if (any(ok)) expression[cbind(ev$i[ok], ev$j[ok])] <- rec$prediction[ok]
  }

  if (nrow(ev)) {
    key <- paste(ev$i, ev$j, sep = ":")
    did <- match(key, paste(det$events$i, det$events$j, sep = ":"))
    events <- det$events[did, , drop = FALSE]
    events$membership <- membership[events$j]
    events$factor_prediction <- rec$factor_prediction
    events$prediction_sd <- rec$prediction_sd
    events$predictability <- rec$predictability
    events$shrinkage <- rec$shrinkage
    events$n_observed_gene <- rec$n_observed_gene
    events$factor_rank <- rec$factor_rank
    events$factor_features <- rec$factor_features
    events$factor_iterations <- rec$factor_iterations
    events$factor_converged <- rec$factor_converged
    events$factor_fold <- if (is.null(rec$factor_fold)) integer(nrow(events)) else rec$factor_fold
    events$bias_calibration <- if (is.null(rec$bias_calibration)) numeric(nrow(events)) else rec$bias_calibration
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
    events$factor_converged <- logical(); events$factor_fold <- integer()
    events$bias_calibration <- numeric(); events$recovery_method <- character()
    events$target_mode <- character(); events$cell_prediction <- numeric()
    events$cell_available <- logical(); events$n_donors <- integer()
    events$bandwidth <- numeric(); events$local_positive_mean <- numeric()
    events$local_positive_variance <- numeric(); events$local_positive_prevalence <- numeric()
    events$effective_donors <- numeric(); events$tree_distance_weighted_mean <- numeric()
    events$embedding_distance_weighted_mean <- numeric(); events$recovered <- numeric()
    events$changed <- logical()
  }

  uncertainty_available <- recovery_method %in% c(
    "p1_stabilized_state", "masked_factor", "tree_local_factor"
  ) && (!any(ok) || all(is.finite(rec$prediction_sd[ok]) & rec$prediction_sd[ok] >= 0))

  predictive_variance <- NULL
  if (uncertainty_available) {
    pv_ok <- ok & is.finite(rec$prediction_sd) & rec$prediction_sd >= 0
    predictive_variance <- .dk_sparse_numeric(
      ev$i[pv_ok], ev$j[pv_ok], rec$prediction_sd[pv_ok]^2,
      nrow(x_work), ncol(x_work), dimnames(x_work)
    )
  }

  detection_scope <- switch(
    detection_method,
    scgacl_gamma_normal = "group",
    alra_global = "all_cells",
    "membership"
  )
  detection_input <- if (identical(detection_method, "scgacl_gamma_normal")) {
    "raw_counts_log1.01"
  } else {
    "recovery_working_matrix"
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
      detection_method = detection_method,
      detection_scope = detection_scope,
      detection_input = detection_input,
      scgacl_dropout_threshold = scgacl_dropout_threshold,
      variance_prior_df = variance_prior_df,
      alra_K = as.integer(alra_K), alra_noise_start = as.integer(alra_noise_start),
      alra_choose_q = as.integer(alra_choose_q), alra_svd_q = as.integer(alra_svd_q),
      recovery_method = recovery_method, factor_target = factor_target,
      factor_rank = factor_rank, factor_features = factor_features,
      factor_ridge = factor_ridge,
      min_feature_observed = min_feature_observed,
      min_target_observed = min_target_observed,
      tree_weight = tree_weight, tree_tau = tree_tau,
      local_k = local_k, candidate_k = candidate_k,
      min_effective_donors = min_effective_donors,
      local_info_kappa = local_info_kappa,
      factor_crossfit_folds = as.integer(factor_crossfit_folds),
      factor_crossfit_seed = as.integer(factor_crossfit_seed),
      support_adaptive_rank = support_adaptive_rank,
      bias_kappa = bias_kappa,
      predictor_smoothing = predictor_smoothing,
      neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
      min_positive_neighbors = min_positive_neighbors,
      neighbor_positive_only = neighbor_positive_only,
      cap_quantile = cap_quantile, seed = seed,
      normalize = normalize,
      normalization_scale_factor = normalization_scale_factor,
      normalization = if (normalize) "ALRA_library_size_log1p" else "none"
    )
  )

  if (return_score) {
    sev <- det$events
    out$score <- .dk_sparse_numeric(
      sev$i, sev$j, sev$confidence,
      det$dimensions[1L], det$dimensions[2L], det$dimnames
    )
    attr(out$score, "zero_only") <- TRUE
    attr(out$score, "detection") <- det
    attr(out$score, "score_type") <- switch(
      detection_method,
      scgacl_gamma_normal = "posterior_dropout_probability",
      alra_global = "binary_native_alra_call",
      "confidence"
    )
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
  if (!is.null(x$settings$detection_scope)) {
    cat(" detection scope:", x$settings$detection_scope, "\n")
  }
  if (!is.null(x$settings$detection_input)) {
    cat(" detection input:", x$settings$detection_input, "\n")
  }
  cat(" recovery engine:", x$settings$recovery_method, "\n")
  if (!is.null(x$settings$normalization)) {
    cat(" normalization:", x$settings$normalization, "\n")
  }
  if (!is.null(x$settings$factor_target) && x$settings$recovery_method %in% c(
    "p1_stabilized_state", "masked_factor", "tree_local_factor"
  )) {
    cat(" recovery target:", x$settings$factor_target, "\n")
  }
  cat(" selected dropout events:", nrow(x$events), "\n")
  cat(" recovered events:", if (nrow(x$events)) sum(x$events$changed) else 0L, "\n")
  invisible(x)
}
