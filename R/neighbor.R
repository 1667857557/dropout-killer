.dk_cell_predict_events <- function(x, embedding, membership, events, k = 30L, sigma = NULL,
                                    min_positive_neighbors = 1L, positive_only = TRUE) {
  n_ev <- nrow(events)
  pred <- rep(NA_real_, n_ev); donors <- integer(n_ev); bandwidth <- rep(NA_real_, n_ev)
  if (!n_ev) return(list(prediction = pred, n_donors = donors, bandwidth = bandwidth))
  k <- max(1L, as.integer(k)); min_positive_neighbors <- max(1L, as.integer(min_positive_neighbors))
  if (!is.null(sigma)) {
    sigma <- as.numeric(sigma)[1L]
    if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive when supplied", call. = FALSE)
  }
  for (m in unique(events$membership)) {
    cells <- which(membership == m); n <- length(cells)
    if (n < 2L) next
    ev_m <- which(events$membership == m)
    query_cells <- sort(unique(events$j[ev_m]))
    query_local <- match(query_cells, cells)
    k_search <- min(n, k + 1L)
    nn <- RANN::nn2(data = embedding[cells, , drop = FALSE],
                    query = embedding[query_cells, , drop = FALSE], k = k_search)
    for (qi in seq_along(query_cells)) {
      cj <- query_cells[qi]; q <- ev_m[events$j[ev_m] == cj]
      ord <- nn$nn.idx[qi, ]; dd <- nn$nn.dists[qi, ]
      keep <- ord != query_local[qi] & is.finite(dd)
      ord <- ord[keep]; dd <- dd[keep]
      if (length(ord) > k) { ord <- ord[seq_len(k)]; dd <- dd[seq_len(k)] }
      if (!length(ord)) next
      sig <- if (is.null(sigma)) {
        posd <- dd[dd > 0]
        if (length(posd)) stats::median(posd) else 1
      } else sigma
      w <- exp(-(dd * dd) / (sig * sig))
      if (!all(is.finite(w)) || sum(w) <= 0) w <- rep(1, length(ord))
      w <- w / sum(w)
      genes <- events$i[q]
      vals <- as.matrix(x[genes, cells[ord], drop = FALSE])
      if (positive_only) {
        obs <- vals > 0
        num <- as.vector(vals %*% w)
        den <- as.vector((obs * 1) %*% w)
        nd <- rowSums(obs)
        p <- num / den
        ok <- is.finite(p) & den > 0 & nd >= min_positive_neighbors
      } else {
        p <- as.vector(vals %*% w); nd <- rep(length(ord), length(q)); ok <- is.finite(p)
      }
      pred[q[ok]] <- p[ok]; donors[q] <- nd; bandwidth[q] <- sig
    }
  }
  list(prediction = pred, n_donors = donors, bandwidth = bandwidth)
}

#' Membership-constrained Gaussian neighbor borrowing
#'
#' @export
weighted_neighbor_prediction <- function(x, embedding, membership, mask, k = 30L, sigma = NULL,
                                         min_positive_neighbors = 1L, positive_only = TRUE,
                                         return_events = FALSE) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  z <- .dk_align_embedding(embedding, nm$cells)
  membership <- .dk_align_membership(membership, nm$cells)
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x)))) stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) events$membership <- membership[events$j] else events$membership <- integer()
  fit <- .dk_cell_predict_events(x, z, membership, events, k, sigma, min_positive_neighbors, positive_only)
  if (return_events) {
    events$prediction <- fit$prediction; events$n_donors <- fit$n_donors; events$bandwidth <- fit$bandwidth
    return(events)
  }
  ok <- is.finite(fit$prediction)
  .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok], nrow(x), ncol(x), list(nm$genes, nm$cells))
}
