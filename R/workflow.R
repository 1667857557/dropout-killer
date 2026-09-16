#' Run the current DropoutKiller workflow
#'
#' Raw RNA counts are library-size normalized to 10,000 and log1p transformed.
#' Recoverable-zero detection is performed only by the calibrated SuperCell
#' hierarchy Lean detector. RNA-only detection rebuilds its geometry from the
#' normalized RNA SVD; paired RNA+ATAC WNN geometry is supplied by the Seurat
#' wrapper. Selected zeros are recovered with P1 stabilized-state recovery by
#' default. Observed nonzero values are never overwritten.
#'
#' @export
dropout_killer <- function(
    x, embedding, group = NULL,
    gamma = 150, k_knn = 5L, rank = "auto",
    neighbor_k = 30L, neighbor_sigma = NULL,
    min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
    cap_quantile = NULL, seed = 12345L, return_score = FALSE,
    recovery_method = c("p1_stabilized_state", "masked_factor",
                        "tree_local_factor", "neighbor"),
    factor_rank = 5L, factor_features = 2000L,
    factor_ridge = 2, min_feature_observed = 20L,
    min_target_observed = 8L,
    detection_method = c("Supercell_hierarchy_Lean_membership",
                         "Supercell_hierarchy_Lean_membership_WNN"),
    factor_target = c("positive", "all_observed"),
    tree_weight = 0.5, tree_tau = NULL,
    local_k = 30L, candidate_k = 100L,
    min_effective_donors = 5,
    local_info_kappa = 5,
    normalize = TRUE, normalization_scale_factor = 1e4,
    factor_crossfit_folds = 5L,
    factor_crossfit_seed = 1L,
    support_adaptive_rank = TRUE,
    bias_kappa = 10,
    predictor_smoothing = 0.25,
    lean_model = NULL, lean_control = list(), lean_geometry = NULL) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  raw_counts <- x
  dimnames(raw_counts) <- list(nm$genes, nm$cells)
  dimnames(x) <- dimnames(raw_counts)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("normalize must be TRUE or FALSE", call. = FALSE)
  if (!is.numeric(normalization_scale_factor) || length(normalization_scale_factor) != 1L ||
      !is.finite(normalization_scale_factor) || normalization_scale_factor <= 0)
    stop("normalization_scale_factor must be a finite value > 0", call. = FALSE)
  if (normalization_scale_factor != 1e4)
    stop("Lean features require normalization_scale_factor=10000", call. = FALSE)
  if (normalize) x <- .dk_alra_library_log(x, scale_factor = normalization_scale_factor)
  z <- .dk_align_embedding(embedding, nm$cells)
  recovery_method <- match.arg(recovery_method)
  detection_method <- match.arg(detection_method)
  factor_target <- match.arg(factor_target)
  group <- .dk_align_vector(group, nm$cells, "group")

  if (is.null(lean_model) && (!is.null(lean_geometry) || !is.null(lean_control$geometry_builder)))
    stop("custom geometry requires a fitted lean_model; calibrate its matching builder explicitly", call. = FALSE)
  if (is.null(lean_geometry)) {
    if (identical(detection_method, .dk_lean_method(TRUE)))
      stop("WNN detection requires build_wnn_supercell geometry or the Seurat wrapper", call. = FALSE)
    detection_embedding <- .dk_lean_svd(x, rank, seed)$embedding
    lean_geometry <- .dk_lean_rna_geometry(detection_embedding, group, gamma, k_knn)
  }
  if (is.null(lean_model)) {
    if (!normalize)
      stop("automatic Lean calibration requires raw counts with normalize=TRUE; otherwise supply lean_model", call. = FALSE)
    if (!is.null(lean_geometry$affinity))
      stop("WNN calibration must rebuild geometry per mask; use the Seurat wrapper or supply lean_model", call. = FALSE)
    lean_model <- do.call(calibrate_lean_detector, c(list(
      counts = raw_counts, group = group, gamma = gamma, k_knn = k_knn,
      neighbor_k = neighbor_k, rank = rank
    ), lean_control))
  }

  membership_fit <- lean_geometry$membership_fit
  membership <- .dk_align_membership(membership_fit$membership, nm$cells)
  membership_fit$membership <- stats::setNames(membership, nm$cells)
  tab <- as.data.frame(table(membership), stringsAsFactors = FALSE)
  names(tab) <- c("membership", "n_cells")
  tab$membership <- as.integer(as.character(tab$membership))
  membership_fit$membership_table <- tab
  lean_geometry$membership_fit <- membership_fit
  hard_recovery_stratum <- membership_fit$cell_stratum[nm$cells]

  det <- supercell_lean_detect(
    x, z, group, model = lean_model, geometry = lean_geometry,
    normalize = FALSE, gamma = gamma, k_knn = k_knn,
    neighbor_k = neighbor_k, rank = rank, seed = seed
  )
  if (!identical(det$settings$detection_method, detection_method))
    stop("requested detector and geometry modality differ", call. = FALSE)
  mask <- .dk_sparse_logical(det$events$i, det$events$j, nrow(x), ncol(x), dimnames(x))
  ev <- .dk_mask_events(mask)
  if (nrow(ev)) ev$membership <- membership[ev$j] else ev$membership <- integer()

  rec <- .dk_recover_events(
    x, ev, membership,
    if (recovery_method %in% c("p1_stabilized_state", "neighbor", "tree_local_factor")) z else NULL,
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
  if (inherits(x, "Matrix")) {
    delta <- .dk_sparse_numeric(ev$i[ok], ev$j[ok], rec$prediction[ok], nrow(x), ncol(x), dimnames(x))
    expression <- x + delta
  } else {
    expression <- x
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

  uncertainty_available <- recovery_method %in% c("p1_stabilized_state", "masked_factor", "tree_local_factor") &&
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
    expression = expression,
    membership = membership,
    membership_fit = membership_fit,
    local_geometry = rec$geometry,
    mask = mask,
    events = events,
    predictive_variance = predictive_variance,
    uncertainty_available = uncertainty_available,
    detection = det,
    settings = list(
      gamma = gamma, k_knn = k_knn, rank = rank,
      detection_method = detection_method,
      detection_scope = "hierarchy_within_broad_group",
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
      sev$i, sev$j, sev$confidence, det$dimensions[1L],
      det$dimensions[2L], det$dimnames
    )
    attr(out$score, "zero_only") <- TRUE
    attr(out$score, "detection") <- det
    attr(out$score, "score_type") <- "balanced_logistic_score_called_events"
  }
  class(out) <- "DropoutKillerResult"
  out
}

#' @export
print.DropoutKillerDetection <- function(x, ...) {
  cat("DropoutKiller detection\n")
  cat(" broad strata:", nrow(x$membership_stats), "\n")
  cat(" called recoverable zeros:", nrow(x$events), "\n")
  cat(" detector:", x$settings$detection_method, "\n")
  invisible(x)
}

#' @export
print.DropoutKillerResult <- function(x, ...) {
  cat("DropoutKiller result\n")
  cat(" dimensions:", paste(dim(x$expression), collapse = " x "), "\n")
  cat(" memberships:", length(unique(x$membership)), "\n")
  cat(" detector:", x$settings$detection_method, "\n")
  cat(" detection scope:", x$settings$detection_scope, "\n")
  cat(" recovery engine:", x$settings$recovery_method, "\n")
  cat(" normalization:", x$settings$normalization, "\n")
  cat(" selected dropout events:", nrow(x$events), "\n")
  cat(" recovered events:", if (nrow(x$events)) sum(x$events$changed) else 0L, "\n")
  invisible(x)
}
