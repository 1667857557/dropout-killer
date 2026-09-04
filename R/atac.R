# ATAC-specific dropout detection and recovery ---------------------------------
# Internal helpers used by the existing Seurat orchestration path.

.dk_atac_binary_matrix <- function(x) {
  if (inherits(x, "Matrix")) {
    y <- x
    if (!length(y@x)) return(Matrix::Matrix(y, sparse = TRUE))
    y@x[] <- 1
    return(Matrix::drop0(y))
  }
  Matrix::Matrix(x > 0, sparse = TRUE)
}

.dk_atac_capture_probability <- function(counts, group = NULL, max_iter = 25L,
                                         tol = 1e-5, eps = 1e-4) {
  counts <- .dk_validate_expression(counts)
  nm <- .dk_names(counts); cells <- nm$cells
  y <- .dk_atac_binary_matrix(counts)
  group <- .dk_align_vector(group, cells, "group")
  if (is.null(group)) group <- rep("all_cells", ncol(counts))
  group <- as.character(group)
  if (anyNA(group) || any(!nzchar(group)))
    stop("ATAC capture group contains missing or empty values", call. = FALSE)
  max_iter <- max(1L, as.integer(max_iter))
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("atac_capture_tol must be > 0", call. = FALSE)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.1)
    stop("atac_capture_eps must be in (0, 0.1)", call. = FALSE)

  q <- rep(NA_real_, ncol(counts)); names(q) <- cells
  lev <- unique(group); diagnostics <- vector("list", length(lev))
  for (kk in seq_along(lev)) {
    g <- lev[kk]; ids <- which(group == g); yg <- y[, ids, drop = FALSE]
    hits <- as.numeric(Matrix::colSums(yg))
    if (!any(hits > 0)) {
      q[ids] <- eps
      diagnostics[[kk]] <- data.frame(group = g, n_cells = length(ids), iterations = 0L,
        converged = TRUE, min_capture = eps, max_capture = eps, stringsAsFactors = FALSE)
      next
    }

    # PIC coordinate descent: q_c = sum_j y_jc / sum_j p_j and
    # p_j = sum_c y_jc / sum_c q_c.  The natural p,q <= 1 boundaries anchor
    # the multiplicative scale without an extra arbitrary depth normalization.
    pcur <- as.numeric(Matrix::rowMeans(yg))
    pcur <- pmin(1 - eps, pmax(eps, pcur))
    qg <- rep(1, length(ids)); converged <- FALSE; last_iter <- 0L
    for (iter in seq_len(max_iter)) {
      denom_p <- sum(pcur)
      if (!is.finite(denom_p) || denom_p <= eps) break
      qnew <- pmin(1 - eps, pmax(eps, hits / denom_p))
      denom_q <- sum(qnew)
      if (!is.finite(denom_q) || denom_q <= eps) break
      pnew <- as.numeric(Matrix::rowSums(yg)) / denom_q
      pnew <- pmin(1 - eps, pmax(eps, pnew))
      last_iter <- iter
      delta <- max(c(abs(qnew - qg), abs(pnew - pcur)))
      qg <- qnew; pcur <- pnew
      if (delta < tol) { converged <- TRUE; break }
    }
    q[ids] <- qg
    diagnostics[[kk]] <- data.frame(group = g, n_cells = length(ids), iterations = last_iter,
      converged = converged, min_capture = min(qg), max_capture = max(qg), stringsAsFactors = FALSE)
  }
  list(capture = q, diagnostics = do.call(rbind, diagnostics), binary = y)
}

.dk_atac_membership_state <- function(binary, membership, capture) {
  cells <- colnames(binary)
  membership <- .dk_align_membership(membership, cells)
  capture <- as.numeric(capture[cells])
  lev <- sort(unique(membership)); map <- match(membership, lev)
  H <- Matrix::sparseMatrix(i = seq_along(map), j = map, x = 1,
    dims = c(length(map), length(lev)), dimnames = list(cells, as.character(lev)))
  n_positive <- binary %*% H
  exposure <- as.numeric(Matrix::crossprod(H, capture))
  exposure2 <- as.numeric(Matrix::crossprod(H, capture * capture))
  names(exposure) <- names(exposure2) <- as.character(lev)
  list(n_positive = n_positive, exposure = exposure, exposure2 = exposure2,
       levels = lev, membership = membership)
}

.dk_atac_bfdr_select <- function(probability, alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("atac_bfdr must be in (0,1)", call. = FALSE)
  probability <- pmin(1, pmax(0, as.numeric(probability)))
  if (!length(probability)) return(list(selected = logical(), bfdr = numeric()))
  ord <- order(probability, decreasing = TRUE)
  q <- cumsum(1 - probability[ord]) / seq_along(ord)
  hit <- which(q <= alpha); k <- if (length(hit)) max(hit) else 0L
  selected <- rep(FALSE, length(probability))
  if (k > 0L) selected[ord[seq_len(k)]] <- TRUE
  bfdr <- rep(NA_real_, length(probability)); bfdr[ord] <- q
  list(selected = selected, bfdr = bfdr)
}

.dk_atac_detect <- function(counts, membership, group = NULL, bfdr = 0.01,
                            pi_min = 0.80, pre_pi = 0.50, score_min = 0.50,
                            min_observed_donors = 2L, min_effective_donors = 5,
                            capture_max_iter = 25L, capture_tol = 1e-5,
                            capture_eps = 1e-4, max_candidates = 5e6) {
  counts <- .dk_validate_expression(counts)
  nm <- .dk_names(counts); cells <- nm$cells; peaks <- nm$genes
  membership <- .dk_align_membership(membership, cells)
  if (!is.numeric(pi_min) || length(pi_min) != 1L || !is.finite(pi_min) || pi_min <= 0 || pi_min >= 1)
    stop("atac_pi_min must be in (0,1)", call. = FALSE)
  if (!is.numeric(pre_pi) || length(pre_pi) != 1L || !is.finite(pre_pi) || pre_pi < 0 || pre_pi >= pi_min)
    stop("atac_pre_pi must be >= 0 and < atac_pi_min", call. = FALSE)
  if (!is.numeric(score_min) || length(score_min) != 1L || !is.finite(score_min) || score_min < 0 || score_min >= 1)
    stop("atac_score_min must be in [0,1)", call. = FALSE)
  min_observed_donors <- max(1L, as.integer(min_observed_donors))
  if (!is.numeric(min_effective_donors) || length(min_effective_donors) != 1L ||
      !is.finite(min_effective_donors) || min_effective_donors < 1)
    stop("atac_min_effective_donors must be >= 1", call. = FALSE)

  cap <- .dk_atac_capture_probability(counts, group = group, max_iter = capture_max_iter,
                                      tol = capture_tol, eps = capture_eps)
  y <- cap$binary; q <- cap$capture
  state <- .dk_atac_membership_state(y, membership, q)
  chunks <- list(); nchunk <- 0L; candidate_total <- 0

  for (mm in seq_along(state$levels)) {
    m <- state$levels[mm]; ids <- which(membership == m)
    if (length(ids) <= 1L) next
    nobs <- as.numeric(state$n_positive[, mm])
    exposure <- state$exposure[mm]; exposure2 <- state$exposure2[mm]
    if (!is.finite(exposure) || exposure <= 0) next
    pi0 <- pmin(1, nobs / exposure)
    eligible <- which(nobs >= min_observed_donors & pi0 >= pre_pi)
    if (!length(eligible)) next

    for (start in seq.int(1L, length(eligible), by = 512L)) {
      pp <- eligible[start:min(length(eligible), start + 511L)]
      b <- as.matrix(y[pp, ids, drop = FALSE]) > 0
      iz <- which(!b, arr.ind = TRUE)
      if (!nrow(iz)) next
      pidx <- pp[iz[, 1L]]; jidx <- ids[iz[, 2L]]; num <- nobs[pidx]
      den <- exposure - q[jidx]; den2 <- exposure2 - q[jidx]^2
      valid <- is.finite(den) & den > capture_eps & is.finite(den2) & den2 > 0
      if (!any(valid)) next
      pi_loo <- rep(0, length(jidx)); neff <- rep(0, length(jidx))
      pi_loo[valid] <- pmin(1, pmax(0, num[valid] / den[valid]))
      neff[valid] <- den[valid]^2 / den2[valid]
      denom <- 1 - pi_loo * q[jidx]; d <- rep(0, length(jidx))
      ok <- valid & denom > .Machine$double.eps
      d[ok] <- pi_loo[ok] * (1 - q[jidx[ok]]) / denom[ok]
      keep <- ok & pi_loo >= pi_min & neff >= min_effective_donors & d >= score_min
      if (!any(keep)) next
      nchunk <- nchunk + 1L
      chunks[[nchunk]] <- data.frame(i = pidx[keep], j = jidx[keep], peak = peaks[pidx[keep]],
        cell = cells[jidx[keep]], membership = m, n_observed_donors = num[keep],
        effective_donors = neff[keep], capture_probability = q[jidx[keep]],
        open_probability = pi_loo[keep], dropout_probability = d[keep], confidence = d[keep],
        stringsAsFactors = FALSE)
      candidate_total <- candidate_total + sum(keep)
      if (is.finite(max_candidates) && candidate_total > max_candidates)
        stop("ATAC candidate count exceeded atac_max_candidates; raise the limit or use stricter atac_pi_min/atac_score_min", call. = FALSE)
    }
  }

  events <- if (length(chunks)) do.call(rbind, chunks) else data.frame(
    i = integer(), j = integer(), peak = character(), cell = character(), membership = integer(),
    n_observed_donors = numeric(), effective_donors = numeric(), capture_probability = numeric(),
    open_probability = numeric(), dropout_probability = numeric(), confidence = numeric(), stringsAsFactors = FALSE)
  sel <- .dk_atac_bfdr_select(events$dropout_probability, bfdr)
  events$bfdr <- sel$bfdr; events$selected <- sel$selected
  chosen <- events[events$selected, , drop = FALSE]
  mask <- .dk_sparse_logical(chosen$i, chosen$j, nrow(counts), ncol(counts), dimnames(counts))
  out <- list(events = events, mask = mask, capture_probability = q,
    capture_diagnostics = cap$diagnostics, state = state, dimensions = dim(counts), dimnames = dimnames(counts),
    settings = list(detection_method = "atac_pic_membership_posterior", bfdr = bfdr,
      pi_min = pi_min, pre_pi = pre_pi, score_min = score_min,
      min_observed_donors = min_observed_donors, min_effective_donors = min_effective_donors,
      capture_max_iter = as.integer(capture_max_iter), capture_tol = capture_tol, capture_eps = capture_eps))
  class(out) <- "DropoutKillerATACDetection"
  out
}

.dk_atac_phi <- function(x, exposure, phi_prior = 1, phi_kappa = 10, phi_floor = 1e-4) {
  out <- rep(phi_prior, nrow(x))
  for (r in seq_len(nrow(x))) {
    pos <- which(x[r, ] > 0 & is.finite(exposure) & exposure > 0)
    n <- length(pos); if (n < 2L) next
    u <- x[r, pos] / exposure[pos]; mu <- mean(u)
    if (!is.finite(mu) || mu <= 0) next
    vv <- stats::var(u); poisson <- mean(mu / exposure[pos])
    ph <- max((vv - poisson) / max(mu^2, 1e-12), phi_floor)
    rr <- n / (n + phi_kappa)
    out[r] <- max(rr * ph + (1 - rr) * phi_prior, phi_floor)
  }
  out
}

.dk_atac_zero_posterior_matrix <- function(binary_block, capture, membership_exposure,
                                           membership_exposure2, capture_eps = 1e-4) {
  nobs <- rowSums(binary_block); den <- membership_exposure - capture
  den2 <- membership_exposure2 - capture^2
  pi <- outer(nobs, den, "/"); pi[!is.finite(pi)] <- 0; pi <- pmin(1, pmax(0, pi))
  d <- pi * matrix(1 - capture, nrow(pi), ncol(pi), byrow = TRUE)
  denom <- 1 - pi * matrix(capture, nrow(pi), ncol(pi), byrow = TRUE)
  ok <- denom > capture_eps; d[ok] <- d[ok] / denom[ok]; d[!ok] <- 0
  s <- d; s[binary_block] <- 1
  list(open = pi, dropout = d, accessibility_weight = s)
}

.dk_atac_recover <- function(counts, detection, membership, embedding,
                             membership_fit = NULL, hard_stratum = NULL,
                             tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
                             candidate_k = 100L, min_effective_donors = 5,
                             local_info_kappa = 5, phi_prior = 1,
                             phi_kappa = 10, phi_floor = 1e-4) {
  counts <- .dk_validate_expression(counts)
  nm <- .dk_names(counts); cells_all <- nm$cells
  membership <- .dk_align_membership(membership, cells_all)
  z <- .dk_align_embedding(embedding, cells_all)
  events <- detection$events[detection$events$selected, , drop = FALSE]
  empty_num <- function() .dk_sparse_numeric(integer(), integer(), numeric(), nrow(counts), ncol(counts), dimnames(counts))
  if (!nrow(events)) return(list(expression = counts, events = events,
    latent_rate = empty_num(), predictive_variance = empty_num(), geometry = NULL))

  geometry <- .dk_build_local_geometry(z, membership, membership_fit, hard_stratum,
    tree_weight = tree_weight, tree_tau = tree_tau, local_k = local_k,
    candidate_k = candidate_k, min_effective_donors = min_effective_donors,
    local_info_kappa = local_info_kappa)
  qcap <- detection$capture_probability[cells_all]
  exposure_cell <- as.numeric(qcap); names(exposure_cell) <- cells_all
  exposure_cell[!is.finite(exposure_cell) | exposure_cell <= 0] <- NA_real_

  E <- nrow(events)
  pred_rate <- pred_count <- pred_sd <- rep(NA_real_, E)
  local_count <- local_exposure <- prior_mean <- prior_phi <- local_neff <- rep(NA_real_, E)
  method <- rep("unavailable", E)

  for (m in unique(events$membership)) {
    qev_m <- which(events$membership == m); ids <- which(membership == m)
    if (length(ids) <= 1L) next
    Wm <- geometry$W[ids, ids, drop = FALSE]
    zero_rows <- which(as.numeric(Matrix::rowSums(Wm)) <= 0)
    if (length(zero_rows)) {
      ii <- rep(zero_rows, each = length(ids) - 1L)
      jj <- unlist(lapply(zero_rows, function(a) setdiff(seq_along(ids), a)), use.names = FALSE)
      Wm <- Wm + Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(length(ids), length(ids)))
    }
    exposure_m <- sum(qcap[ids]); exposure2_m <- sum(qcap[ids]^2)
    p_all <- unique(events$i[qev_m])

    for (start in seq.int(1L, length(p_all), by = 256L)) {
      pp <- p_all[start:min(length(p_all), start + 255L)]
      X <- as.matrix(counts[pp, ids, drop = FALSE]); B <- X > 0
      post0 <- .dk_atac_zero_posterior_matrix(B, qcap[ids], exposure_m, exposure2_m,
                                               detection$settings$capture_eps)
      S <- post0$accessibility_weight; ex <- exposure_cell[ids]
      invalid <- !is.finite(ex) | ex <= 0
      if (any(invalid)) { S[, invalid] <- 0; ex[invalid] <- 0 }
      weighted_count <- as.matrix(X %*% Matrix::t(Wm))
      weighted_exposure <- as.matrix(sweep(S, 2L, ex, "*") %*% Matrix::t(Wm))
      W2m <- Wm * Wm
      den <- as.matrix(S %*% Matrix::t(Wm)); den2 <- as.matrix((S * S) %*% Matrix::t(W2m))
      neff <- matrix(0, nrow(S), ncol(S)); okn <- den2 > 0; neff[okn] <- den[okn]^2 / den2[okn]
      total_exp <- rowSums(sweep(S, 2L, ex, "*")); total_count <- rowSums(X)
      phi <- .dk_atac_phi(X, ex, phi_prior = phi_prior, phi_kappa = phi_kappa, phi_floor = phi_floor)

      qev <- qev_m[events$i[qev_m] %in% pp]
      rr <- match(events$i[qev], pp); cc <- match(events$j[qev], ids); ij <- cbind(rr, cc)
      wc <- weighted_count[ij]; we <- weighted_exposure[ij]
      target_exp <- S[ij] * exposure_cell[events$j[qev]]
      mu <- total_count[rr] / pmax(total_exp[rr] - target_exp, .Machine$double.eps)
      ph <- phi[rr]; local_n <- neff[ij]
      available <- is.finite(wc) & wc >= 0 & is.finite(we) & we > 0 & is.finite(mu) & mu > 0 &
        is.finite(ph) & ph > 0 & is.finite(local_n) & local_n >= min_effective_donors &
        is.finite(exposure_cell[events$j[qev]]) & exposure_cell[events$j[qev]] > 0
      if (!any(available)) next
      qa <- qev[available]
      aa <- 1 / ph[available] + wc[available]
      bb <- 1 / (ph[available] * mu[available]) + we[available]
      theta <- aa / bb; vv <- aa / (bb * bb)
      pred_rate[qa] <- theta
      pred_count[qa] <- exposure_cell[events$j[qa]] * theta
      pred_sd[qa] <- exposure_cell[events$j[qa]] * sqrt(pmax(vv, 0))
      local_count[qa] <- wc[available]; local_exposure[qa] <- we[available]
      prior_mean[qa] <- mu[available]; prior_phi[qa] <- ph[available]
      local_neff[qa] <- local_n[available]; method[qa] <- "atac_local_gamma_poisson_eb"
    }
  }

  changed <- is.finite(pred_count) & pred_count > 0
  if (inherits(counts, "Matrix")) {
    expression <- counts + .dk_sparse_numeric(events$i[changed], events$j[changed], pred_count[changed],
      nrow(counts), ncol(counts), dimnames(counts))
  } else {
    expression <- counts
    if (any(changed)) expression[cbind(events$i[changed], events$j[changed])] <- pred_count[changed]
  }
  latent_rate <- .dk_sparse_numeric(events$i[changed], events$j[changed], pred_rate[changed],
    nrow(counts), ncol(counts), dimnames(counts))
  predictive_variance <- .dk_sparse_numeric(events$i[changed], events$j[changed], pred_sd[changed]^2,
    nrow(counts), ncol(counts), dimnames(counts))
  events$local_weighted_count <- local_count; events$local_accessible_exposure <- local_exposure
  events$prior_mean <- prior_mean; events$prior_phi <- prior_phi
  events$recovered_rate <- pred_rate; events$recovered <- pred_count; events$prediction_sd <- pred_sd
  events$effective_recovery_donors <- local_neff; events$recovery_method <- method; events$changed <- changed
  list(expression = expression, events = events, latent_rate = latent_rate,
       predictive_variance = predictive_variance, geometry = geometry, capture_exposure = exposure_cell)
}

.dk_atac_dropout_killer <- function(counts, embedding, membership, membership_fit = NULL,
    group = NULL, hard_stratum = NULL, bfdr = 0.01, pi_min = 0.80, pre_pi = 0.50,
    score_min = 0.50, min_observed_donors = 2L, min_effective_donors = 5,
    capture_max_iter = 25L, capture_tol = 1e-5, capture_eps = 1e-4, max_candidates = 5e6,
    tree_weight = 0.5, tree_tau = NULL, local_k = 30L, candidate_k = 100L,
    local_info_kappa = 5, phi_prior = 1, phi_kappa = 10, phi_floor = 1e-4) {
  det <- .dk_atac_detect(counts, membership = membership, group = group, bfdr = bfdr,
    pi_min = pi_min, pre_pi = pre_pi, score_min = score_min,
    min_observed_donors = min_observed_donors, min_effective_donors = min_effective_donors,
    capture_max_iter = capture_max_iter, capture_tol = capture_tol, capture_eps = capture_eps,
    max_candidates = max_candidates)
  rec <- .dk_atac_recover(counts, det, membership = membership, embedding = embedding,
    membership_fit = membership_fit, hard_stratum = hard_stratum, tree_weight = tree_weight,
    tree_tau = tree_tau, local_k = local_k, candidate_k = candidate_k,
    min_effective_donors = min_effective_donors, local_info_kappa = local_info_kappa,
    phi_prior = phi_prior, phi_kappa = phi_kappa, phi_floor = phi_floor)
  out <- list(expression = rec$expression, latent_rate = rec$latent_rate, mask = det$mask,
    events = rec$events, candidate_events = det$events, predictive_variance = rec$predictive_variance,
    capture_probability = det$capture_probability, capture_diagnostics = det$capture_diagnostics,
    detection = det, local_geometry = rec$geometry,
    settings = c(det$settings, list(recovery_method = "atac_local_gamma_poisson_eb",
      tree_weight = tree_weight, tree_tau = tree_tau, local_k = local_k, candidate_k = candidate_k,
      local_info_kappa = local_info_kappa, phi_prior = phi_prior, phi_kappa = phi_kappa,
      phi_floor = phi_floor)))
  class(out) <- "DropoutKillerATACResult"
  out
}
