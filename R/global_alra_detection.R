.dk_alra_rsvd <- function(A, k, q) {
  k <- as.integer(k)
  q <- as.integer(q)
  if (k < 1L || k > min(dim(A))) {
    stop("ALRA SVD rank must be in [1, min(dim(A))]", call. = FALSE)
  }
  if (q < 0L) stop("ALRA power iterations must be >= 0", call. = FALSE)
  rsvd::rsvd(A, k = k, q = q)
}

.dk_alra_choose_k_native <- function(A, K = 100L, thresh = 6,
                                     noise_start = 80L, q = 2L,
                                     seed = NULL) {
  # Direct mathematical port of KlugerLab/ALRA::choose_k(). Unlike the former
  # DropoutKiller cell-class implementation, K and noise_start are NOT adapted
  # to small blocks: the original ALRA validity checks are retained verbatim.
  K <- as.integer(K)
  noise_start <- as.integer(noise_start)
  if (K > min(dim(A))) {
    stop("For an m by n matrix, K must not exceed min(m,n), matching original ALRA.",
         call. = FALSE)
  }
  if (noise_start > K - 5L) {
    stop("There need to be at least 5 singular values considered noise, matching original ALRA.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  noise_svals <- noise_start:K
  fit <- .dk_alra_rsvd(A, K, q = q)
  d <- fit$d
  diffs <- d[seq_len(length(d) - 1L)] - d[2:length(d)]
  mu <- mean(diffs[noise_svals - 1L])
  sigma <- stats::sd(diffs[noise_svals - 1L])
  num_of_sds <- (diffs - mu) / sigma
  hit <- which(num_of_sds > thresh)
  if (!length(hit)) {
    stop("Original ALRA choose_k found no singular-value spacing above the threshold.",
         call. = FALSE)
  }
  k <- max(hit)
  list(k = k, num_of_sds = num_of_sds, d = d,
       K = K, noise_start = noise_start)
}

.dk_alra_native_block <- function(y, ids = seq_len(ncol(y)),
                                  block_id = 1L, block_label = "all_cells",
                                  rank = "auto", quantile_prob = 0.001,
                                  K = 100L, rank_z = 6, noise_start = 80L,
                                  choose_q = 2L, svd_q = 10L, seed = 1L) {
  # ALRA stores cells in rows and genes in columns.
  A <- t(as.matrix(y))
  mind <- min(dim(A))
  if (mind < 2L) return(NULL)

  numeric_rank <- is.numeric(rank) && length(rank) == 1L &&
    is.finite(rank) && rank >= 1
  if (is.numeric(rank) && !numeric_rank) {
    stop("numeric rank must be a finite value >= 1", call. = FALSE)
  }

  # The released ALRA code relies on one R RNG stream: choose_k() consumes the
  # first randomized SVD and alra() then immediately performs the final q=10
  # randomized SVD. Set the seed once here rather than reseeding each stage.
  set.seed(as.integer(seed))
  if (numeric_rank) {
    k <- as.integer(rank)
    if (k > mind) stop("numeric ALRA rank exceeds min(matrix dimensions)", call. = FALSE)
    choice <- list(k = k, K = k, noise_start = NA_integer_, d = numeric())
  } else {
    if (!identical(rank, "auto")) {
      stop("rank must be a positive integer or 'auto'", call. = FALSE)
    }
    choice <- .dk_alra_choose_k_native(
      A, K = K, thresh = rank_z, noise_start = noise_start,
      q = choose_q, seed = NULL
    )
    k <- choice$k
  }

  fit <- .dk_alra_rsvd(A, k = k, q = svd_q)
  Uscaled <- sweep(
    fit$u[, seq_len(k), drop = FALSE], 2L,
    fit$d[seq_len(k)], FUN = "*"
  )
  lr <- Uscaled %*% t(fit$v[, seq_len(k), drop = FALSE])

  # Exact threshold used by the released ALRA code:
  # abs(apply(A_norm_rank_k, 2, quantile, quantile.prob)). R's default
  # quantile type (type=7) is intentionally retained.
  tau <- abs(apply(
    lr, 2L, stats::quantile,
    probs = quantile_prob, names = FALSE
  ))

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

.dk_original_alra_detect <- function(x, rank = "auto",
                                     quantile_prob = 0.001,
                                     seed = 1L, K = 100L, rank_z = 6,
                                     noise_start = 80L, choose_q = 2L,
                                     svd_q = 10L) {
  x <- .dk_validate_expression(x)
  nm <- .dk_names(x)
  if (!is.numeric(quantile_prob) || length(quantile_prob) != 1L ||
      quantile_prob <= 0 || quantile_prob >= 0.5) {
    stop("quantile_prob must be between 0 and 0.5", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2) {
    stop("alra_K must be a finite value >= 2", call. = FALSE)
  }
  if (!is.numeric(rank_z) || length(rank_z) != 1L ||
      !is.finite(rank_z) || rank_z <= 0) {
    stop("rank_z must be a finite value > 0", call. = FALSE)
  }
  if (!is.numeric(noise_start) || length(noise_start) != 1L ||
      !is.finite(noise_start) || noise_start < 2) {
    stop("alra_noise_start must be a finite value >= 2", call. = FALSE)
  }
  if (!is.numeric(choose_q) || length(choose_q) != 1L ||
      !is.finite(choose_q) || choose_q < 0) {
    stop("alra_choose_q must be a finite value >= 0", call. = FALSE)
  }
  if (!is.numeric(svd_q) || length(svd_q) != 1L ||
      !is.finite(svd_q) || svd_q < 0) {
    stop("alra_svd_q must be a finite value >= 0", call. = FALSE)
  }

  fit <- .dk_alra_native_block(
    x, ids = seq_len(ncol(x)), block_id = 1L, block_label = "all_cells",
    rank = rank, quantile_prob = quantile_prob, K = K,
    rank_z = rank_z, noise_start = noise_start,
    choose_q = choose_q, svd_q = svd_q, seed = seed
  )
  if (is.null(fit)) stop("ALRA requires at least two cells and two genes", call. = FALSE)

  ev <- if (!is.null(fit$events) && nrow(fit$events)) {
    fit$events
  } else {
    data.frame(
      i = integer(), j = integer(), gene = character(), cell = character(),
      membership = integer(), detection_block = character(), lowrank = numeric(),
      threshold = numeric(), null_sigma = numeric(), z_score = numeric(),
      p_value = numeric(), q_value = numeric(), confidence = numeric(),
      confidence_fallback = logical(), variance_weight = numeric(),
      alra_margin = numeric(), stringsAsFactors = FALSE
    )
  }

  st <- data.frame(
    membership = 1L,
    detection_block = "all_cells",
    n_cells = ncol(x),
    rank = fit$rank,
    status = "ok",
    alra_candidates = fit$n_candidates,
    zero_tests = fit$n_zero,
    prior_sigma = NA_real_,
    detection_method = "alra_global",
    K = fit$K,
    noise_start = fit$noise_start,
    stringsAsFactors = FALSE
  )

  out <- list(
    events = ev,
    membership_stats = st,
    dimensions = dim(x),
    dimnames = list(nm$genes, nm$cells),
    settings = list(
      rank = rank,
      quantile_prob = quantile_prob,
      seed = seed,
      detection_method = "alra_global",
      K = as.integer(K),
      rank_z = rank_z,
      noise_start = as.integer(noise_start),
      choose_q = as.integer(choose_q),
      svd_q = as.integer(svd_q),
      detection_scope = "all_cells"
    )
  )
  class(out) <- "DropoutKillerDetection"
  out
}

# Internal backwards-compatible symbol. The historical function name implied
# group-wise blocks; from v0.8 it intentionally delegates to the original ALRA
# all-cell implementation so no code path labelled ALRA silently uses a modified
# cell-class algorithm.
.dk_global_alra_detect <- function(x, group = NULL, rank = "auto",
                                   quantile_prob = 0.001, min_cells = 8L,
                                   seed = 1L, K = 100L, rank_z = 6,
                                   noise_start = 80L, choose_q = 2L,
                                   svd_q = 10L) {
  .dk_original_alra_detect(
    x, rank = rank, quantile_prob = quantile_prob,
    seed = seed, K = K, rank_z = rank_z,
    noise_start = noise_start, choose_q = choose_q, svd_q = svd_q
  )
}
