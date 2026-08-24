.dk_svd_fit <- function(y, rank = "auto", max_rank = 20L, rank_z = 6, seed = 1L) {
  nr <- nrow(y); nc <- ncol(y); mind <- min(nr, nc)
  if (mind < 2L) return(NULL)
  if (!is.numeric(max_rank) || length(max_rank) != 1L || !is.finite(max_rank) || max_rank < 1) stop("max_rank must be a finite value >= 1", call. = FALSE)
  if (!is.numeric(rank_z) || length(rank_z) != 1L || !is.finite(rank_z) || rank_z <= 0) stop("rank_z must be a finite value > 0", call. = FALSE)
  max_k <- max(1L, min(as.integer(max_rank), mind - 1L))
  numeric_rank <- is.numeric(rank) && length(rank) == 1L && is.finite(rank) && rank >= 1
  if (is.numeric(rank) && !numeric_rank) stop("numeric rank must be a finite value >= 1", call. = FALSE)
  K <- if (numeric_rank) max(1L, min(as.integer(rank), mind - 1L)) else {
    if (!identical(rank, "auto")) stop("rank must be a positive integer or 'auto'", call. = FALSE)
    max_k
  }
  if (nc <= nr && nc <= 256L) {
    gram <- as.matrix(crossprod(y))
    eig <- eigen(gram, symmetric = TRUE)
    d_all <- sqrt(pmax(eig$values, 0))
    V_all <- eig$vectors
    fit_type <- "gram_exact"
  } else {
    set.seed(seed)
    fit <- irlba::irlba(y, nu = K, nv = K)
    d_all <- fit$d; V_all <- fit$v
    fit_type <- "irlba"
  }
  d <- d_all[seq_len(min(K, length(d_all)))]
  if (numeric_rank) {
    k <- min(K, length(d))
  } else {
    diffs <- if (length(d) > 1L) d[-length(d)] - d[-1L] else numeric()
    if (length(diffs) < 6L) {
      k <- 1L
    } else {
      noise_n <- min(20L, max(5L, floor(length(diffs) / 3L)))
      tail_d <- utils::tail(diffs, noise_n)
      s <- stats::sd(tail_d); mu <- mean(tail_d)
      if (!is.finite(s) || s <= 0) k <- 1L else {
        z <- (diffs - mu) / s
        hit <- which(is.finite(z) & z > rank_z)
        k <- if (length(hit)) max(hit) else 1L
      }
    }
    k <- max(1L, min(as.integer(k), K, ncol(V_all)))
  }
  V <- V_all[, seq_len(k), drop = FALSE]
  lowrank <- as.matrix((y %*% V) %*% t(V))
  list(lowrank = lowrank, rank = k, singular_values = d_all, fit_type = fit_type)
}

#' Detect ALRA-gated zero events within memberships
#'
#' The adaptive threshold follows ALRA's negative-tail idea. A zero is eligible
#' only when its low-rank reconstruction exceeds the magnitude of the lower
#' gene-specific quantile. A one-sided confidence against a symmetric Gaussian
#' biological-zero error model is then estimated from negative reconstructed
#' values. The confidence is evidence against the biological-zero null, not a
#' calibrated Bayesian posterior probability of dropout.
#'
#' @export
local_alra_detect <- function(x, membership, rank = "auto", max_rank = 20L, rank_z = 6,
                              quantile_prob = 0.001, min_cells = 8L, min_negative = 3L,
                              seed = 1L) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  if (!is.numeric(quantile_prob) || length(quantile_prob) != 1L || quantile_prob <= 0 || quantile_prob >= 0.5)
    stop("quantile_prob must be between 0 and 0.5", call. = FALSE)
  min_cells <- max(2L, as.integer(min_cells)); min_negative <- max(1L, as.integer(min_negative))
  events <- list(); stats_out <- list(); e <- 0L
  lev <- sort(unique(membership))
  for (ii in seq_along(lev)) {
    m <- lev[ii]; ids <- which(membership == m); n <- length(ids)
    if (n < min_cells) {
      stats_out[[ii]] <- data.frame(membership = m, n_cells = n, rank = NA_integer_, status = "too_small", alra_candidates = 0L)
      next
    }
    y <- x[, ids, drop = FALSE]
    fit <- .dk_svd_fit(y, rank = rank, max_rank = max_rank, rank_z = rank_z, seed = seed + ii - 1L)
    if (is.null(fit)) {
      stats_out[[ii]] <- data.frame(membership = m, n_cells = n, rank = NA_integer_, status = "svd_unavailable", alra_candidates = 0L)
      next
    }
    yh <- fit$lowrank
    qv <- matrixStats::rowQuantiles(yh, probs = quantile_prob, drop = TRUE)
    tau <- abs(qv)
    supported <- is.finite(qv) & qv < 0
    tau[!supported] <- Inf
    zero <- as.matrix(y == 0)
    pass <- zero & sweep(yh, 1L, tau, FUN = ">")
    idx <- which(pass, arr.ind = TRUE)
    if (!nrow(idx)) {
      stats_out[[ii]] <- data.frame(membership = m, n_cells = n, rank = fit$rank, status = "ok", alra_candidates = 0L)
      next
    }
    neg <- yh < 0
    nneg <- rowSums(neg)
    negsq <- yh; negsq[!neg] <- 0
    sigma <- sqrt(rowSums(negsq * negsq) / pmax(nneg, 1L))
    rr <- idx[, 1L]; cc <- idx[, 2L]; val <- yh[idx]
    conf <- stats::pnorm(val / sigma[rr])
    fallback <- nneg[rr] < min_negative | !is.finite(sigma[rr]) | sigma[rr] <= 0 | !is.finite(conf)
    conf[fallback] <- 0.5
    conf <- pmax(0, pmin(1, conf))
    e <- e + 1L
    events[[e]] <- data.frame(i = rr, j = ids[cc], membership = m, lowrank = val,
                              threshold = tau[rr], confidence = conf,
                              confidence_fallback = fallback, stringsAsFactors = FALSE)
    stats_out[[ii]] <- data.frame(membership = m, n_cells = n, rank = fit$rank, status = "ok", alra_candidates = nrow(idx))
  }
  ev <- if (length(events)) do.call(rbind, events) else data.frame(i = integer(), j = integer(), membership = integer(), lowrank = numeric(), threshold = numeric(), confidence = numeric(), confidence_fallback = logical())
  if (nrow(ev)) {
    ev$gene <- nm$genes[ev$i]; ev$cell <- nm$cells[ev$j]
    ev <- ev[, c("i", "j", "gene", "cell", "membership", "lowrank", "threshold", "confidence", "confidence_fallback")]
  }
  st <- if (length(stats_out)) do.call(rbind, stats_out) else data.frame()
  out <- list(events = ev, membership_stats = st, dimensions = dim(x), dimnames = list(nm$genes, nm$cells),
              settings = list(rank = rank, max_rank = max_rank, rank_z = rank_z, quantile_prob = quantile_prob,
                              min_cells = min_cells, min_negative = min_negative, seed = seed))
  class(out) <- "DropoutKillerDetection"
  out
}

#' Sparse local ALRA confidence matrix
#'
#' @export
local_alra_score <- function(x, membership, rank = "auto", max_rank = 20L, rank_z = 6,
                             quantile_prob = 0.001, min_cells = 8L, min_negative = 3L, seed = 1L) {
  det <- local_alra_detect(x, membership, rank, max_rank, rank_z, quantile_prob, min_cells, min_negative, seed)
  ev <- det$events
  s <- .dk_sparse_numeric(ev$i, ev$j, ev$confidence, det$dimensions[1L], det$dimensions[2L], det$dimnames)
  attr(s, "zero_only") <- TRUE
  attr(s, "detection") <- det
  s
}

#' Select high-confidence dropout mask
#'
#' @export
select_dropout_mask <- function(score, threshold = 0.95, x = NULL) {
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold <= 0 || threshold > 1) stop("threshold must be in (0,1]", call. = FALSE)
  if (inherits(score, "DropoutKillerDetection")) {
    ev <- score$events; keep <- ev$confidence >= threshold
    return(.dk_sparse_logical(ev$i[keep], ev$j[keep], score$dimensions[1L], score$dimensions[2L], score$dimnames))
  }
  if (inherits(score, "sparseMatrix")) {
    sm <- Matrix::summary(score); keep <- sm$x >= threshold
    i <- sm$i[keep]; j <- sm$j[keep]
    nr <- nrow(score); nc <- ncol(score); dn <- dimnames(score)
  } else {
    sc <- as.matrix(score); idx <- which(sc >= threshold, arr.ind = TRUE)
    i <- if (nrow(idx)) idx[, 1L] else integer(); j <- if (nrow(idx)) idx[, 2L] else integer()
    nr <- nrow(sc); nc <- ncol(sc); dn <- dimnames(sc)
  }
  if (!is.null(x) && length(i)) {
    x <- .dk_validate_expression(x)
    if (!identical(dim(x), c(nr, nc))) stop("x and score dimensions differ", call. = FALSE)
    keep <- as.vector(x[cbind(i, j)] == 0); i <- i[keep]; j <- j[keep]
  }
  .dk_sparse_logical(i, j, nr, nc, dn)
}

#' @export
print.DropoutKillerDetection <- function(x, ...) {
  cat("DropoutKiller detection\n")
  cat(" memberships:", nrow(x$membership_stats), "\n")
  cat(" ALRA-gated zero events:", nrow(x$events), "\n")
  if (nrow(x$membership_stats)) cat(" skipped small memberships:", sum(x$membership_stats$status == "too_small"), "\n")
  invisible(x)
}
