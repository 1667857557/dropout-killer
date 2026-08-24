.dk_recover_events <- function(x, events, membership, embedding,
                               neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only,
                               cap_quantile = NULL) {
  if (!is.logical(neighbor_positive_only) || length(neighbor_positive_only) != 1L || is.na(neighbor_positive_only)) stop("neighbor_positive_only must be TRUE or FALSE", call. = FALSE)
  cp <- .dk_cell_predict_events(x, embedding, membership, events, neighbor_k, neighbor_sigma,
                                min_positive_neighbors, positive_only = neighbor_positive_only)
  available <- is.finite(cp$prediction)
  pred <- cp$prediction
  pred[!available | pred < 0] <- 0
  if (!is.null(cap_quantile)) {
    if (!is.numeric(cap_quantile) || length(cap_quantile) != 1L || !is.finite(cap_quantile) || cap_quantile <= 0.5 || cap_quantile > 1) stop("cap_quantile must be in (0.5,1]", call. = FALSE)
    for (m in unique(events$membership)) {
      q <- which(events$membership == m); cells <- which(membership == m)
      for (g in unique(events$i[q])) {
        z <- q[events$i[q] == g]; obs <- as.numeric(x[g, cells]); obs <- obs[obs > 0]
        if (length(obs)) pred[z] <- pmin(pred[z], as.numeric(stats::quantile(obs, cap_quantile, names = FALSE)))
      }
    }
  }
  list(prediction = pred, cell_prediction = cp$prediction, cell_available = available,
       n_donors = cp$n_donors, bandwidth = cp$bandwidth)
}

#' Selectively recover masked zero events
#'
#' Recovery uses only membership-constrained Gaussian neighbor borrowing. The
#' donor weights are estimated from latent-space distances and renormalized over
#' eligible donors for each target gene.
#'
#' @export
recover_dropout_expression <- function(x, mask, membership, embedding,
                                       neighbor_k = 30L, neighbor_sigma = NULL,
                                       min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                                       cap_quantile = NULL, return_details = FALSE) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  z <- .dk_align_embedding(embedding, nm$cells)
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x)))) stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) {
    if (any(events$i > nrow(x) | events$j > ncol(x))) stop("mask dimensions exceed x", call. = FALSE)
    is_zero <- as.vector(x[cbind(events$i, events$j)] == 0)
    if (!all(is_zero)) stop("mask contains observed non-zero entries; selective recovery refuses to overwrite them", call. = FALSE)
    events$membership <- membership[events$j]
  } else events$membership <- integer()
  fit <- .dk_recover_events(x, events, membership, z, neighbor_k, neighbor_sigma,
                            min_positive_neighbors, neighbor_positive_only, cap_quantile)
  ok <- fit$prediction > 0 & is.finite(fit$prediction)
  if (inherits(x, "Matrix")) {
    delta <- .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok], nrow(x), ncol(x), dimnames(x))
    out <- x + delta
  } else {
    out <- x
    if (any(ok)) out[cbind(events$i[ok], events$j[ok])] <- fit$prediction[ok]
  }
  events$cell_prediction <- fit$cell_prediction
  events$cell_available <- fit$cell_available
  events$n_donors <- fit$n_donors
  events$bandwidth <- fit$bandwidth
  events$recovered <- fit$prediction
  events$changed <- ok
  if (return_details) list(expression = out, events = events) else out
}

#' Run selective dropout detection and recovery
#'
#' High-confidence zero events are detected from membership-local low-rank
#' structure and recovered only from membership-local cells using Gaussian
#' latent-space weights estimated from cell distances.
#'
#' @export
dropout_killer <- function(x, embedding, membership = NULL, group = NULL, split_by = NULL,
                           gamma = 20, k_knn = 5L, approximate = "auto", approx_n = 20000L,
                           rank = "auto", max_rank = 20L, rank_z = 6, quantile_prob = 0.001,
                           threshold = 0.95, min_cells = 8L, min_negative = 3L,
                           neighbor_k = 30L, neighbor_sigma = NULL,
                           min_positive_neighbors = 1L, neighbor_positive_only = TRUE,
                           cap_quantile = NULL, seed = 12345L, return_score = FALSE) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  z <- .dk_align_embedding(embedding, nm$cells)
  group <- .dk_align_vector(group, nm$cells, "group")
  split_by <- .dk_align_vector(split_by, nm$cells, "split_by")
  membership_fit <- NULL
  if (is.null(membership)) {
    membership_fit <- build_supercell_membership(z, group = group, split_by = split_by, gamma = gamma,
                                                  k_knn = k_knn, approximate = approximate,
                                                  approx_n = approx_n, seed = seed)
    membership <- membership_fit$membership
  } else membership <- .dk_align_membership(membership, nm$cells)
  det <- local_alra_detect(x, membership, rank = rank, max_rank = max_rank, rank_z = rank_z,
                           quantile_prob = quantile_prob, min_cells = min_cells,
                           min_negative = min_negative, seed = seed)
  mask <- select_dropout_mask(det, threshold = threshold)
  ev <- .dk_mask_events(mask)
  if (nrow(ev)) ev$membership <- membership[ev$j] else ev$membership <- integer()
  rec <- .dk_recover_events(x, ev, membership, z, neighbor_k, neighbor_sigma,
                            min_positive_neighbors, neighbor_positive_only, cap_quantile)
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
    events$cell_prediction <- rec$cell_prediction
    events$cell_available <- rec$cell_available
    events$n_donors <- rec$n_donors
    events$bandwidth <- rec$bandwidth
    events$recovered <- rec$prediction
    events$changed <- ok
  } else {
    events <- det$events[FALSE, , drop = FALSE]
    events$cell_prediction <- numeric(); events$cell_available <- logical()
    events$n_donors <- integer(); events$bandwidth <- numeric()
    events$recovered <- numeric(); events$changed <- logical()
  }
  out <- list(expression = expression, membership = membership, membership_fit = membership_fit,
              mask = mask, events = events, detection = det,
              settings = list(gamma = gamma, k_knn = k_knn, approximate = approximate, approx_n = approx_n,
                              rank = rank, max_rank = max_rank, rank_z = rank_z, quantile_prob = quantile_prob,
                              threshold = threshold, min_cells = min_cells, min_negative = min_negative,
                              neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
                              min_positive_neighbors = min_positive_neighbors,
                              neighbor_positive_only = neighbor_positive_only,
                              cap_quantile = cap_quantile, seed = seed))
  if (return_score) {
    sev <- det$events
    out$score <- .dk_sparse_numeric(sev$i, sev$j, sev$confidence, det$dimensions[1L], det$dimensions[2L], det$dimnames)
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
  cat(" high-confidence dropout events:", nrow(x$events), "\n")
  cat(" recovered events:", if (nrow(x$events)) sum(x$events$changed) else 0L, "\n")
  invisible(x)
}
