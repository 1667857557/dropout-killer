.dk_recover_events <- function(x, events, membership, embedding, alpha, adjacency,
                               neighbor_k, neighbor_sigma, min_positive_neighbors, neighbor_positive_only,
                               gene_standardize, cap_quantile = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) stop("alpha must be in [0,1]", call. = FALSE)
  if (!is.logical(neighbor_positive_only) || length(neighbor_positive_only) != 1L || is.na(neighbor_positive_only)) stop("neighbor_positive_only must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(gene_standardize) || length(gene_standardize) != 1L || is.na(gene_standardize)) stop("gene_standardize must be TRUE or FALSE", call. = FALSE)
  cp <- .dk_cell_predict_events(x, embedding, membership, events, neighbor_k, neighbor_sigma,
                                min_positive_neighbors, positive_only = neighbor_positive_only)
  gp <- .dk_gene_predict_events(x, adjacency, membership, events, standardize = gene_standardize)
  c_ok <- is.finite(cp$prediction); g_ok <- is.finite(gp$prediction)
  cw <- alpha * c_ok; gw <- (1 - alpha) * g_ok; den <- cw + gw
  pred <- rep(0, nrow(events)); use <- den > 0
  pred[use] <- (ifelse(c_ok[use], alpha * cp$prediction[use], 0) + ifelse(g_ok[use], (1 - alpha) * gp$prediction[use], 0)) / den[use]
  pred[!is.finite(pred) | pred < 0] <- 0
  if (!is.null(cap_quantile)) {
    if (!is.numeric(cap_quantile) || cap_quantile <= 0.5 || cap_quantile > 1) stop("cap_quantile must be in (0.5,1]", call. = FALSE)
    for (m in unique(events$membership)) {
      q <- which(events$membership == m); cells <- which(membership == m)
      for (g in unique(events$i[q])) {
        z <- q[events$i[q] == g]; obs <- as.numeric(x[g, cells]); obs <- obs[obs > 0]
        if (length(obs)) pred[z] <- pmin(pred[z], as.numeric(stats::quantile(obs, cap_quantile, names = FALSE)))
      }
    }
  }
  list(prediction = pred, cell_prediction = cp$prediction, gene_prediction = gp$prediction,
       cell_available = c_ok, gene_available = g_ok, n_donors = cp$n_donors,
       bandwidth = cp$bandwidth, network_support = gp$support)
}

#' Selectively recover masked zero events
#'
#' @export
recover_dropout_expression <- function(x, mask, membership, embedding, alpha = 0.75,
                                       adjacency = NULL, neighbor_k = 30L, neighbor_sigma = NULL,
                                       min_positive_neighbors = 1L, neighbor_positive_only = TRUE, gene_standardize = TRUE,
                                       cap_quantile = NULL, return_details = FALSE) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  z <- .dk_align_embedding(embedding, nm$cells)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) stop("alpha must be in [0,1]", call. = FALSE)
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x)))) stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) {
    if (any(events$i > nrow(x) | events$j > ncol(x))) stop("mask dimensions exceed x", call. = FALSE)
    is_zero <- as.vector(x[cbind(events$i, events$j)] == 0)
    if (!all(is_zero)) stop("mask contains observed non-zero entries; selective recovery refuses to overwrite them", call. = FALSE)
    events$membership <- membership[events$j]
  } else events$membership <- integer()
  A <- if (is.null(adjacency)) NULL else prepare_gene_network(adjacency, nm$genes)
  fit <- .dk_recover_events(x, events, membership, z, alpha, A, neighbor_k, neighbor_sigma,
                            min_positive_neighbors, neighbor_positive_only, gene_standardize, cap_quantile)
  ok <- fit$prediction > 0 & is.finite(fit$prediction)
  if (inherits(x, "Matrix")) {
    delta <- .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok], nrow(x), ncol(x), dimnames(x))
    out <- x + delta
  } else {
    out <- x
    if (any(ok)) out[cbind(events$i[ok], events$j[ok])] <- fit$prediction[ok]
  }
  events$cell_prediction <- fit$cell_prediction
  events$gene_prediction <- fit$gene_prediction
  events$cell_available <- fit$cell_available
  events$gene_available <- fit$gene_available
  events$n_donors <- fit$n_donors
  events$bandwidth <- fit$bandwidth
  events$network_support <- fit$network_support
  events$recovered <- fit$prediction
  events$changed <- ok
  if (return_details) list(expression = out, events = events) else out
}

#' Run selective dropout detection and recovery
#'
#' @export
dropout_killer <- function(x, embedding, membership = NULL, group = NULL, split_by = NULL,
                           gamma = 20, k_knn = 5L, approximate = "auto", approx_n = 20000L,
                           rank = "auto", max_rank = 20L, rank_z = 6, quantile_prob = 0.001,
                           threshold = 0.95, min_cells = 8L, min_negative = 3L,
                           alpha = 0.75, adjacency = NULL, gene_networks = NULL,
                           gene_network_weights = NULL, neighbor_k = 30L, neighbor_sigma = NULL,
                           min_positive_neighbors = 1L, neighbor_positive_only = TRUE, gene_standardize = TRUE,
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
  if (!is.null(adjacency) && !is.null(gene_networks)) stop("supply adjacency or gene_networks, not both", call. = FALSE)
  A <- if (!is.null(gene_networks)) combine_gene_prior(gene_networks, nm$genes, gene_network_weights) else if (!is.null(adjacency)) prepare_gene_network(adjacency, nm$genes) else NULL
  det <- local_alra_detect(x, membership, rank = rank, max_rank = max_rank, rank_z = rank_z,
                           quantile_prob = quantile_prob, min_cells = min_cells,
                           min_negative = min_negative, seed = seed)
  mask <- select_dropout_mask(det, threshold = threshold)
  ev <- .dk_mask_events(mask)
  if (nrow(ev)) ev$membership <- membership[ev$j] else ev$membership <- integer()
  rec <- .dk_recover_events(x, ev, membership, z, alpha, A, neighbor_k, neighbor_sigma,
                            min_positive_neighbors, neighbor_positive_only, gene_standardize, cap_quantile)
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
    events$gene_prediction <- rec$gene_prediction
    events$cell_available <- rec$cell_available
    events$gene_available <- rec$gene_available
    events$n_donors <- rec$n_donors
    events$bandwidth <- rec$bandwidth
    events$network_support <- rec$network_support
    events$recovered <- rec$prediction
    events$changed <- ok
  } else {
    events <- det$events[FALSE, , drop = FALSE]
    events$cell_prediction <- numeric(); events$gene_prediction <- numeric(); events$cell_available <- logical(); events$gene_available <- logical()
    events$n_donors <- integer(); events$bandwidth <- numeric(); events$network_support <- integer(); events$recovered <- numeric(); events$changed <- logical()
  }
  out <- list(expression = expression, membership = membership, membership_fit = membership_fit,
              mask = mask, events = events, detection = det,
              settings = list(gamma = gamma, k_knn = k_knn, approximate = approximate, approx_n = approx_n,
                              rank = rank, max_rank = max_rank, rank_z = rank_z, quantile_prob = quantile_prob,
                              threshold = threshold, min_cells = min_cells, min_negative = min_negative,
                              alpha = alpha, neighbor_k = neighbor_k, neighbor_sigma = neighbor_sigma,
                              min_positive_neighbors = min_positive_neighbors, neighbor_positive_only = neighbor_positive_only,
                              gene_standardize = gene_standardize,
                              cap_quantile = cap_quantile, gene_prior = !is.null(A), seed = seed))
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
