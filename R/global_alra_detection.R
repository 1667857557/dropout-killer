.dk_alra_rsvd <- function(A, k, q, seed) {
  k <- as.integer(k); q <- as.integer(q); seed <- as.integer(seed)
  if (k < 1L || k >= min(dim(A))) stop("ALRA SVD rank must be in [1, min(dim(A))-1]", call. = FALSE)
  if (q < 0L) stop("ALRA power iterations must be >= 0", call. = FALSE)
  set.seed(seed)
  rsvd::rsvd(A, k = k, q = q)
}

.dk_alra_choose_k_native <- function(A, K = 100L, thresh = 6,
                                     noise_start = 80L, q = 2L, seed = 1L) {
  mind <- min(dim(A))
  if (mind < 3L) return(list(k = 1L, d = numeric(), K = 1L,
                             noise_start = NA_integer_, num_of_sds = numeric()))
  K <- min(as.integer(K), mind - 1L)
  if (K < 2L) return(list(k = 1L, d = numeric(), K = K,
                          noise_start = NA_integer_, num_of_sds = numeric()))
  fit <- .dk_alra_rsvd(A, K, q = q, seed = seed)
  d <- fit$d
  diffs <- d[-length(d)] - d[-1L]
  if (length(diffs) < 5L) {
    return(list(k = 1L, d = d, K = K, noise_start = NA_integer_,
                num_of_sds = rep(NA_real_, length(diffs))))
  }
  noise_start_eff <- min(as.integer(noise_start), K - 4L)
  noise_start_eff <- max(2L, noise_start_eff)
  noise_idx <- (noise_start_eff:K) - 1L
  noise_idx <- noise_idx[noise_idx >= 1L & noise_idx <= length(diffs)]
  mu <- mean(diffs[noise_idx])
  sigma <- stats::sd(diffs[noise_idx])
  if (!is.finite(sigma) || sigma <= 0) {
    z <- rep(NA_real_, length(diffs)); k <- 1L
  } else {
    z <- (diffs - mu) / sigma
    hit <- which(is.finite(z) & z > thresh)
    k <- if (length(hit)) max(hit) else 1L
  }
  k <- max(1L, min(as.integer(k), K))
  list(k = k, d = d, K = K, noise_start = noise_start_eff, num_of_sds = z)
}

.dk_alra_native_block <- function(y, ids, block_id, block_label,
                                  rank = "auto", quantile_prob = 0.001,
                                  K = 100L, rank_z = 6, noise_start = 80L,
                                  choose_q = 2L, svd_q = 10L, seed = 1L) {
  A <- t(as.matrix(y))
  mind <- min(dim(A))
  if (mind < 2L) return(NULL)
  numeric_rank <- is.numeric(rank) && length(rank) == 1L && is.finite(rank) && rank >= 1
  if (is.numeric(rank) && !numeric_rank) stop("numeric rank must be a finite value >= 1", call. = FALSE)
  if (numeric_rank) {
    k <- max(1L, min(as.integer(rank), mind - 1L))
    choice <- list(k = k, K = k, noise_start = NA_integer_, d = numeric())
  } else {
    if (!identical(rank, "auto")) stop("rank must be a positive integer or 'auto'", call. = FALSE)
    choice <- .dk_alra_choose_k_native(
      A, K = K, thresh = rank_z, noise_start = noise_start,
      q = choose_q, seed = seed
    )
    k <- choice$k
  }
  fit <- .dk_alra_rsvd(A, k = k, q = svd_q, seed = seed + 100000L)
  Uscaled <- sweep(fit$u[, seq_len(k), drop = FALSE], 2L,
                   fit$d[seq_len(k)], FUN = "*")
  lr <- Uscaled %*% t(fit$v[, seq_len(k), drop = FALSE])
  tau <- abs(apply(lr, 2L, stats::quantile, probs = quantile_prob,
                   names = FALSE))
  zero <- A == 0
  pass <- zero & sweep(lr, 2L, tau, FUN = ">")
  idx <- which(pass, arr.ind = TRUE)
  if (!nrow(idx)) {
    ev <- NULL
  } else {
    cc_local <- idx[, 1L]
    gg <- idx[, 2L]
    val <- lr[idx]
    margin <- val - tau[gg]
    ev <- data.frame(
      i = gg,
      j = ids[cc_local],
      gene = rownames(y)[gg],
      cell = colnames(y)[cc_local],
      membership = NA_integer_,
      detection_block = block_label,
      lowrank = val,
      threshold = tau[gg],
      null_sigma = NA_real_,
      z_score = NA_real_,
      p_value = NA_real_,
      q_value = NA_real_,
      confidence = 1,
      confidence_fallback = FALSE,
      variance_weight = NA_real_,
      alra_margin = margin,
      stringsAsFactors = FALSE
    )
  }
  list(
    events = ev,
    rank = k,
    K = choice$K,
    noise_start = choice$noise_start,
    n_zero = sum(zero),
    n_candidates = if (is.null(ev)) 0L else nrow(ev),
    block_id = block_id,
    block_label = block_label
  )
}

.dk_global_alra_detect <- function(x, group = NULL, rank = "auto",
                                   quantile_prob = 0.001, min_cells = 8L,
                                   seed = 1L, K = 100L, rank_z = 6,
                                   noise_start = 80L, choose_q = 2L,
                                   svd_q = 10L) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  if (!is.numeric(quantile_prob) || length(quantile_prob) != 1L ||
      quantile_prob <= 0 || quantile_prob >= 0.5)
    stop("quantile_prob must be between 0 and 0.5", call. = FALSE)
  min_cells <- max(2L, as.integer(min_cells))
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2)
    stop("alra_K must be a finite value >= 2", call. = FALSE)
  if (!is.numeric(rank_z) || length(rank_z) != 1L || !is.finite(rank_z) || rank_z <= 0)
    stop("rank_z must be a finite value > 0", call. = FALSE)
  if (!is.numeric(noise_start) || length(noise_start) != 1L || !is.finite(noise_start) || noise_start < 2)
    stop("alra_noise_start must be a finite value >= 2", call. = FALSE)
  if (!is.numeric(choose_q) || length(choose_q) != 1L || !is.finite(choose_q) || choose_q < 0)
    stop("alra_choose_q must be a finite value >= 0", call. = FALSE)
  if (!is.numeric(svd_q) || length(svd_q) != 1L || !is.finite(svd_q) || svd_q < 0)
    stop("alra_svd_q must be a finite value >= 0", call. = FALSE)
  if (is.null(group)) {
    grp <- rep("all_cells", ncol(x))
  } else {
    grp <- as.character(group)
    if (length(grp) != ncol(x)) stop("group length must equal number of cells", call. = FALSE)
    if (!is.null(names(group)) && !is.null(nm$cells)) {
      hit <- match(nm$cells, names(group))
      if (anyNA(hit)) stop("named group does not cover all cells", call. = FALSE)
      grp <- as.character(group[hit])
    }
    if (anyNA(grp) || any(!nzchar(grp))) stop("group contains missing or empty cell classes", call. = FALSE)
  }
  lev <- unique(grp)
  events <- list(); stats_out <- vector("list", length(lev)); e <- 0L
  for (ii in seq_along(lev)) {
    label <- lev[ii]
    ids <- which(grp == label)
    n <- length(ids)
    if (n < min_cells) {
      stats_out[[ii]] <- data.frame(
        membership = ii, detection_block = label, n_cells = n,
        rank = NA_integer_, status = "too_small", alra_candidates = 0L,
        zero_tests = 0L, prior_sigma = NA_real_,
        detection_method = "alra_global_by_group", K = NA_integer_,
        noise_start = NA_integer_, stringsAsFactors = FALSE
      )
      next
    }
    y <- x[, ids, drop = FALSE]
    fit <- .dk_alra_native_block(
      y, ids = ids, block_id = ii, block_label = label,
      rank = rank, quantile_prob = quantile_prob, K = K,
      rank_z = rank_z, noise_start = noise_start,
      choose_q = choose_q, svd_q = svd_q, seed = seed + ii - 1L
    )
    if (is.null(fit)) {
      stats_out[[ii]] <- data.frame(
        membership = ii, detection_block = label, n_cells = n,
        rank = NA_integer_, status = "svd_unavailable", alra_candidates = 0L,
        zero_tests = 0L, prior_sigma = NA_real_,
        detection_method = "alra_global_by_group", K = NA_integer_,
        noise_start = NA_integer_, stringsAsFactors = FALSE
      )
      next
    }
    if (!is.null(fit$events) && nrow(fit$events)) {
      e <- e + 1L; events[[e]] <- fit$events
    }
    stats_out[[ii]] <- data.frame(
      membership = ii, detection_block = label, n_cells = n,
      rank = fit$rank, status = "ok", alra_candidates = fit$n_candidates,
      zero_tests = fit$n_zero, prior_sigma = NA_real_,
      detection_method = "alra_global_by_group", K = fit$K,
      noise_start = fit$noise_start, stringsAsFactors = FALSE
    )
  }
  ev <- if (length(events)) do.call(rbind, events) else data.frame(
    i = integer(), j = integer(), gene = character(), cell = character(),
    membership = integer(), detection_block = character(), lowrank = numeric(),
    threshold = numeric(), null_sigma = numeric(), z_score = numeric(),
    p_value = numeric(), q_value = numeric(), confidence = numeric(),
    confidence_fallback = logical(), variance_weight = numeric(),
    alra_margin = numeric(), stringsAsFactors = FALSE
  )
  st <- if (length(stats_out)) do.call(rbind, stats_out) else data.frame()
  out <- list(
    events = ev, membership_stats = st, dimensions = dim(x),
    dimnames = list(nm$genes, nm$cells),
    settings = list(
      rank = rank, quantile_prob = quantile_prob, min_cells = min_cells,
      seed = seed, detection_method = "alra_global_by_group",
      K = as.integer(K), rank_z = rank_z,
      noise_start = as.integer(noise_start),
      choose_q = as.integer(choose_q), svd_q = as.integer(svd_q),
      detection_scope = if (is.null(group)) "all_cells" else "group"
    )
  )
  class(out) <- "DropoutKillerDetection"
  out
}
