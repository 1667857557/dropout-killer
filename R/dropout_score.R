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

.dk_shrunk_negative_null <- function(yh, min_negative = 3L, variance_prior_df = 10) {
  min_negative <- max(1L, as.integer(min_negative))
  if (!is.numeric(variance_prior_df) || length(variance_prior_df) != 1L ||
      !is.finite(variance_prior_df) || variance_prior_df < 0)
    stop("variance_prior_df must be a finite value >= 0", call. = FALSE)
  neg <- yh < 0
  nneg <- rowSums(neg)
  negsq <- yh
  negsq[!neg] <- 0
  ss <- rowSums(negsq * negsq)
  local_sigma2 <- ss / pmax(nneg, 1L)
  valid_local <- nneg >= min_negative & is.finite(local_sigma2) & local_sigma2 > 0
  prior_sigma2 <- if (any(valid_local)) stats::median(local_sigma2[valid_local], na.rm = TRUE) else NA_real_
  if (!is.finite(prior_sigma2) || prior_sigma2 <= 0) {
    total_n <- sum(nneg)
    prior_sigma2 <- if (total_n > 0L) sum(ss) / total_n else NA_real_
  }
  if (!is.finite(prior_sigma2) || prior_sigma2 <= 0) {
    return(list(sigma = rep(NA_real_, nrow(yh)), sigma2 = rep(NA_real_, nrow(yh)),
                n_negative = nneg, prior_sigma = NA_real_, shrinkage_weight = numeric(nrow(yh)),
                supported = rep(FALSE, nrow(yh))))
  }
  local_sigma2[!is.finite(local_sigma2) | local_sigma2 <= 0] <- prior_sigma2
  w <- nneg / (nneg + variance_prior_df)
  if (variance_prior_df == 0) w[] <- as.numeric(nneg > 0L)
  sigma2 <- w * local_sigma2 + (1 - w) * prior_sigma2
  supported <- is.finite(sigma2) & sigma2 > 0
  list(sigma = sqrt(pmax(sigma2, 0)), sigma2 = sigma2, n_negative = nneg,
       prior_sigma = sqrt(prior_sigma2), shrinkage_weight = w, supported = supported)
}

.dk_gene_bh <- function(p_value, gene_index, n_zero_by_gene) {
  q <- rep(1, length(p_value))
  if (!length(p_value)) return(q)
  by_gene <- split(seq_along(p_value), gene_index)
  for (ix in by_gene) {
    g <- gene_index[ix[1L]]
    q[ix] <- stats::p.adjust(p_value[ix], method = "BH", n = max(length(ix), n_zero_by_gene[g]))
  }
  pmax(0, pmin(1, q))
}

.dk_detect_eb_zero_null <- function(y, yh, ids, membership_id, fit, min_negative, variance_prior_df) {
  null <- .dk_shrunk_negative_null(yh, min_negative = min_negative,
                                   variance_prior_df = variance_prior_df)
  zero <- as.matrix(y == 0)
  n_zero_by_gene <- rowSums(zero)
  idx <- which(zero & yh > 0, arr.ind = TRUE)
  if (!nrow(idx)) {
    return(list(events = NULL, n_candidates = 0L, n_zero = sum(n_zero_by_gene),
                prior_sigma = null$prior_sigma))
  }
  rr <- idx[, 1L]; cc <- idx[, 2L]
  val <- yh[idx]
  sigma <- null$sigma[rr]
  z_score <- val / sigma
  p_value <- stats::pnorm(z_score, lower.tail = FALSE)
  bad <- !is.finite(z_score) | !is.finite(p_value) | !null$supported[rr]
  p_value[bad] <- 1
  q_value <- .dk_gene_bh(p_value, rr, n_zero_by_gene)
  confidence <- 1 - q_value

  # Sparse detection output contains only positions with more evidence for the
  # expressed/dropout alternative than for the symmetric zero null. Selection
  # at the default 0.95 threshold is therefore gene-wise BH q <= 0.05.
  keep <- is.finite(confidence) & confidence > 0.5
  if (!any(keep)) {
    return(list(events = NULL, n_candidates = 0L, n_zero = sum(n_zero_by_gene),
                prior_sigma = null$prior_sigma))
  }
  rr <- rr[keep]; cc <- cc[keep]; val <- val[keep]; sigma <- sigma[keep]
  z_score <- z_score[keep]; p_value <- p_value[keep]; q_value <- q_value[keep]
  confidence <- confidence[keep]
  nneg <- null$n_negative[rr]
  w <- null$shrinkage_weight[rr]
  fallback <- nneg < min_negative
  ev <- data.frame(
    i = rr, j = ids[cc], membership = membership_id, lowrank = val,
    threshold = NA_real_, null_sigma = sigma, z_score = z_score,
    p_value = p_value, q_value = q_value, confidence = confidence,
    confidence_fallback = fallback, variance_weight = w,
    stringsAsFactors = FALSE
  )
  list(events = ev, n_candidates = nrow(ev), n_zero = sum(n_zero_by_gene),
       prior_sigma = null$prior_sigma)
}

.dk_detect_alra_quantile <- function(y, yh, ids, membership_id, fit,
                                     quantile_prob, min_negative) {
  qv <- matrixStats::rowQuantiles(yh, probs = quantile_prob, drop = TRUE)
  tau <- abs(qv)
  supported <- is.finite(qv) & qv < 0
  tau[!supported] <- Inf
  zero <- as.matrix(y == 0)
  pass <- zero & sweep(yh, 1L, tau, FUN = ">")
  idx <- which(pass, arr.ind = TRUE)
  if (!nrow(idx))
    return(list(events = NULL, n_candidates = 0L, n_zero = sum(zero), prior_sigma = NA_real_))
  neg <- yh < 0
  nneg <- rowSums(neg)
  negsq <- yh; negsq[!neg] <- 0
  sigma <- sqrt(rowSums(negsq * negsq) / pmax(nneg, 1L))
  rr <- idx[, 1L]; cc <- idx[, 2L]; val <- yh[idx]
  conf <- stats::pnorm(val / sigma[rr])
  fallback <- nneg[rr] < min_negative | !is.finite(sigma[rr]) | sigma[rr] <= 0 | !is.finite(conf)
  conf[fallback] <- 0.5
  conf <- pmax(0, pmin(1, conf))
  ev <- data.frame(
    i = rr, j = ids[cc], membership = membership_id, lowrank = val,
    threshold = tau[rr], null_sigma = sigma[rr], z_score = val / sigma[rr],
    p_value = 1 - conf, q_value = NA_real_, confidence = conf,
    confidence_fallback = fallback, variance_weight = NA_real_,
    stringsAsFactors = FALSE
  )
  list(events = ev, n_candidates = nrow(ev), n_zero = sum(zero), prior_sigma = NA_real_)
}

#' Detect low-rank-supported zero events within memberships
#'
#' The default finite-sample detector keeps ALRA's low-rank reconstruction but
#' replaces the empirical 0.1% tail gate with a symmetric biological-zero null.
#' Gene-specific zero-null variance is estimated from negative reconstructed
#' values and shrunk toward a robust membership-level variance center. Positive
#' reconstructed values at observed zeros receive one-sided Gaussian p-values;
#' Benjamini-Hochberg adjustment is performed separately for each gene within
#' each membership. `confidence` is `1 - q_value`, not a Bayesian posterior
#' probability. Thus the default `threshold = 0.95` corresponds to gene-wise
#' BH-adjusted q <= 0.05 under the working null model.
#'
#' The historical ALRA empirical-quantile gate remains available with
#' `detection_method = "alra_quantile"` for reproducibility.
#'
#' @export
local_alra_detect <- function(x, membership, rank = "auto", max_rank = 20L, rank_z = 6,
                              quantile_prob = 0.001, min_cells = 8L, min_negative = 3L,
                              seed = 1L,
                              detection_method = c("eb_zero_null", "alra_quantile"),
                              variance_prior_df = 10) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  detection_method <- match.arg(detection_method)
  if (!is.numeric(quantile_prob) || length(quantile_prob) != 1L || quantile_prob <= 0 || quantile_prob >= 0.5)
    stop("quantile_prob must be between 0 and 0.5", call. = FALSE)
  min_cells <- max(2L, as.integer(min_cells)); min_negative <- max(1L, as.integer(min_negative))
  if (!is.numeric(variance_prior_df) || length(variance_prior_df) != 1L ||
      !is.finite(variance_prior_df) || variance_prior_df < 0)
    stop("variance_prior_df must be a finite value >= 0", call. = FALSE)
  events <- list(); stats_out <- list(); e <- 0L
  lev <- sort(unique(membership))
  for (ii in seq_along(lev)) {
    m <- lev[ii]; ids <- which(membership == m); n <- length(ids)
    if (n < min_cells) {
      stats_out[[ii]] <- data.frame(
        membership = m, n_cells = n, rank = NA_integer_, status = "too_small",
        alra_candidates = 0L, zero_tests = 0L, prior_sigma = NA_real_,
        detection_method = detection_method
      )
      next
    }
    y <- x[, ids, drop = FALSE]
    fit <- .dk_svd_fit(y, rank = rank, max_rank = max_rank, rank_z = rank_z,
                       seed = seed + ii - 1L)
    if (is.null(fit)) {
      stats_out[[ii]] <- data.frame(
        membership = m, n_cells = n, rank = NA_integer_, status = "svd_unavailable",
        alra_candidates = 0L, zero_tests = 0L, prior_sigma = NA_real_,
        detection_method = detection_method
      )
      next
    }
    yh <- fit$lowrank
    detm <- if (detection_method == "eb_zero_null") {
      .dk_detect_eb_zero_null(y, yh, ids, m, fit, min_negative, variance_prior_df)
    } else {
      .dk_detect_alra_quantile(y, yh, ids, m, fit, quantile_prob, min_negative)
    }
    if (!is.null(detm$events) && nrow(detm$events)) {
      e <- e + 1L
      events[[e]] <- detm$events
    }
    stats_out[[ii]] <- data.frame(
      membership = m, n_cells = n, rank = fit$rank, status = "ok",
      alra_candidates = detm$n_candidates, zero_tests = detm$n_zero,
      prior_sigma = detm$prior_sigma, detection_method = detection_method
    )
  }
  ev <- if (length(events)) do.call(rbind, events) else data.frame(
    i = integer(), j = integer(), membership = integer(), lowrank = numeric(),
    threshold = numeric(), null_sigma = numeric(), z_score = numeric(),
    p_value = numeric(), q_value = numeric(), confidence = numeric(),
    confidence_fallback = logical(), variance_weight = numeric()
  )
  if (nrow(ev)) {
    ev$gene <- nm$genes[ev$i]; ev$cell <- nm$cells[ev$j]
    ev <- ev[, c("i", "j", "gene", "cell", "membership", "lowrank", "threshold",
                 "null_sigma", "z_score", "p_value", "q_value", "confidence",
                 "confidence_fallback", "variance_weight")]
  }
  st <- if (length(stats_out)) do.call(rbind, stats_out) else data.frame()
  out <- list(
    events = ev, membership_stats = st, dimensions = dim(x),
    dimnames = list(nm$genes, nm$cells),
    settings = list(
      rank = rank, max_rank = max_rank, rank_z = rank_z,
      quantile_prob = quantile_prob, min_cells = min_cells,
      min_negative = min_negative, seed = seed,
      detection_method = detection_method,
      variance_prior_df = variance_prior_df
    )
  )
  class(out) <- "DropoutKillerDetection"
  out
}

#' Sparse local low-rank zero-null confidence matrix
#'
#' @export
local_alra_score <- function(x, membership, rank = "auto", max_rank = 20L, rank_z = 6,
                             quantile_prob = 0.001, min_cells = 8L, min_negative = 3L,
                             seed = 1L,
                             detection_method = c("eb_zero_null", "alra_quantile"),
                             variance_prior_df = 10) {
  det <- local_alra_detect(
    x, membership, rank, max_rank, rank_z, quantile_prob, min_cells,
    min_negative, seed, detection_method, variance_prior_df
  )
  ev <- det$events
  s <- .dk_sparse_numeric(ev$i, ev$j, ev$confidence, det$dimensions[1L],
                          det$dimensions[2L], det$dimnames)
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
  cat(" low-rank-supported zero events:", nrow(x$events), "\n")
  cat(" detector:", if (!is.null(x$settings$detection_method)) x$settings$detection_method else "legacy", "\n")
  if (nrow(x$membership_stats)) cat(" skipped small memberships:", sum(x$membership_stats$status == "too_small"), "\n")
  invisible(x)
}
