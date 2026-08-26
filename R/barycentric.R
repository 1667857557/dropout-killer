.dk_project_simplex <- function(v) {
  v <- as.numeric(v)
  n <- length(v)
  if (!n) return(v)
  if (n == 1L) return(1)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  rho <- which(u - cssv / seq_along(u) > 0)
  if (!length(rho)) return(rep(1 / n, n))
  rho <- max(rho)
  theta <- cssv[rho] / rho
  w <- pmax(v - theta, 0)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) rep(1 / n, n) else w / s
}

.dk_barycentric_control <- function(membership_penalty = 1, barycentric_lambda = 1,
                                    barycentric_iter = 20L, barycentric_tol = 1e-6,
                                    local_k = 30L, candidate_k = 100L,
                                    min_effective_donors = 5) {
  if (!is.numeric(membership_penalty) || length(membership_penalty) != 1L ||
      !is.finite(membership_penalty) || membership_penalty < 0)
    stop("membership_penalty must be a non-negative finite number", call. = FALSE)
  if (!is.numeric(barycentric_lambda) || length(barycentric_lambda) != 1L ||
      !is.finite(barycentric_lambda) || barycentric_lambda <= 0)
    stop("barycentric_lambda must be a positive finite number", call. = FALSE)
  if (!is.numeric(barycentric_iter) || length(barycentric_iter) != 1L ||
      !is.finite(barycentric_iter) || barycentric_iter < 1)
    stop("barycentric_iter must be >= 1", call. = FALSE)
  if (!is.numeric(barycentric_tol) || length(barycentric_tol) != 1L ||
      !is.finite(barycentric_tol) || barycentric_tol <= 0)
    stop("barycentric_tol must be a positive finite number", call. = FALSE)
  for (z in list(local_k = local_k, candidate_k = candidate_k, min_effective_donors = min_effective_donors)) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z < 1)
      stop("local neighborhood controls must be >= 1", call. = FALSE)
  }
  list(membership_penalty = as.numeric(membership_penalty),
       barycentric_lambda = as.numeric(barycentric_lambda),
       barycentric_iter = as.integer(barycentric_iter),
       barycentric_tol = as.numeric(barycentric_tol),
       local_k = as.integer(local_k), candidate_k = as.integer(candidate_k),
       min_effective_donors = as.numeric(min_effective_donors))
}

.dk_soft_recovery_stratum <- function(cells, membership_fit = NULL, hard_stratum = NULL) {
  if (!is.null(hard_stratum)) {
    s <- .dk_align_vector(hard_stratum, cells, "hard_stratum")
    return(factor(as.character(s), levels = unique(as.character(s))))
  }
  has_hard <- !is.null(membership_fit) && isTRUE(membership_fit$settings$has_hard_stratum)
  if (has_hard && !is.null(membership_fit$cell_stratum)) {
    s <- membership_fit$cell_stratum
    if (!is.null(names(s))) s <- s[cells]
    if (length(s) == length(cells) && !anyNA(s))
      return(factor(as.character(s), levels = unique(as.character(s))))
  }
  factor(rep("all", length(cells)), levels = "all")
}

.dk_membership_tree_pair_distance <- function(tree_index, cells, membership,
                                              query_membership, donor_membership) {
  n <- length(query_membership)
  out <- rep(NA_real_, n)
  if (!n || is.null(tree_index)) return(out)
  same <- query_membership == donor_membership
  out[same] <- 0
  use <- which(!same & !is.na(query_membership) & !is.na(donor_membership))
  if (!length(use)) return(out)
  available <- names(tree_index$ancestors)
  if (!length(available)) return(out)
  by_membership <- split(cells[cells %in% available], membership[cells %in% available])
  representative <- vapply(by_membership, function(v) v[1L], character(1))
  a <- pmin(query_membership[use], donor_membership[use])
  b <- pmax(query_membership[use], donor_membership[use])
  key <- paste(a, b, sep = ":")
  ukey <- unique(key)
  parts <- strsplit(ukey, ":", fixed = TRUE)
  ua <- vapply(parts, `[`, character(1), 1L)
  ub <- vapply(parts, `[`, character(1), 2L)
  ra <- unname(representative[ua]); rb <- unname(representative[ub])
  good <- !is.na(ra) & !is.na(rb)
  d <- rep(NA_real_, length(ukey))
  if (any(good)) d[good] <- .dk_tree_distance(tree_index, ra[good], rb[good])
  names(d) <- ukey
  out[use] <- unname(d[key])
  out
}

.dk_barycentric_weights <- function(delta, prior, lambda = 1, max_iter = 20L, tol = 1e-6) {
  k <- nrow(delta)
  if (!k) return(list(weight = numeric(), iterations = 0L, converged = FALSE,
                      state_error = NA_real_, prior_error = NA_real_))
  if (k == 1L) {
    e <- sum(delta[1L, ]^2)
    return(list(weight = 1, iterations = 0L, converged = TRUE,
                state_error = e, prior_error = e))
  }
  prior <- .dk_project_simplex(prior)
  w <- prior
  step <- 1 / (2 * (sum(delta * delta) + lambda))
  converged <- FALSE
  iterations <- 0L
  for (it in seq_len(max_iter)) {
    residual <- as.vector(crossprod(w, delta))
    grad <- 2 * as.vector(delta %*% residual) + 2 * lambda * (w - prior)
    w_new <- .dk_project_simplex(w - step * grad)
    iterations <- it
    if (max(abs(w_new - w)) <= tol) {
      w <- w_new
      converged <- TRUE
      break
    }
    w <- w_new
  }
  w[w < 1e-10] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) w <- prior else w <- w / sw
  state_error <- sum(as.vector(crossprod(w, delta))^2)
  prior_error <- sum(as.vector(crossprod(prior, delta))^2)
  list(weight = w, iterations = iterations, converged = converged,
       state_error = state_error, prior_error = prior_error)
}

.dk_build_barycentric_geometry <- function(embedding, membership, query_cells,
                                           membership_fit = NULL, hard_stratum = NULL,
                                           tree_weight = 0.5, tree_tau = NULL,
                                           local_k = 30L, candidate_k = 100L,
                                           membership_penalty = 1,
                                           barycentric_lambda = 1,
                                           barycentric_iter = 20L,
                                           barycentric_tol = 1e-6,
                                           min_effective_donors = 5) {
  ctl <- .dk_barycentric_control(membership_penalty, barycentric_lambda,
                                 barycentric_iter, barycentric_tol,
                                 local_k, candidate_k, min_effective_donors)
  if (!is.numeric(tree_weight) || length(tree_weight) != 1L || !is.finite(tree_weight) ||
      tree_weight < 0 || tree_weight > 1)
    stop("tree_weight must be in [0,1]", call. = FALSE)
  if (!is.null(tree_tau) && (!is.numeric(tree_tau) || length(tree_tau) != 1L ||
      !is.finite(tree_tau) || tree_tau <= 0))
    stop("tree_tau must be NULL or a positive finite number", call. = FALSE)
  z <- as.matrix(embedding)
  n <- nrow(z); cells <- rownames(z)
  if (is.null(cells)) stop("embedding must have cell row names for barycentric recovery", call. = FALSE)
  query_cells <- unique(as.integer(query_cells))
  query_cells <- query_cells[is.finite(query_cells) & query_cells >= 1L & query_cells <= n]
  nq <- length(query_cells)
  strata <- .dk_soft_recovery_stratum(cells, membership_fit, hard_stratum)
  if (!nq) {
    W <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(0L, n),
                              dimnames = list(character(), cells))
    return(list(W = W, W2 = W, prior = W, prior2 = W, query_cells = integer(),
                bandwidth = numeric(), candidate_count = integer(), state_error = numeric(),
                prior_error = numeric(), tree_distance_weighted_mean = numeric(),
                embedding_distance_weighted_mean = numeric(), same_membership_weight = numeric(),
                iterations = integer(), converged = logical(), cell_stratum = as.character(strata),
                control = ctl))
  }
  qrow <- setNames(seq_len(nq), query_cells)
  bandwidth <- rep(NA_real_, nq); candidate_count <- integer(nq)
  state_error <- rep(NA_real_, nq); prior_error <- rep(NA_real_, nq)
  tree_mean <- rep(NA_real_, nq); embed_mean <- rep(NA_real_, nq)
  same_membership_weight <- rep(NA_real_, nq); iterations <- integer(nq); converged <- logical(nq)
  wi <- integer(); wj <- integer(); wx <- numeric(); px <- numeric()
  qstrata <- as.character(strata[query_cells])
  for (s in unique(qstrata)) {
    qg <- query_cells[qstrata == s]
    data_idx <- which(as.character(strata) == s)
    ns <- length(data_idx)
    if (ns <= 1L) next
    k_use <- min(ns - 1L, ctl$candidate_k)
    nn <- RANN::nn2(data = z[data_idx, , drop = FALSE], query = z[qg, , drop = FALSE],
                    k = min(ns, k_use + 1L))
    cand <- vector("list", length(qg)); dist <- vector("list", length(qg)); hvec <- numeric(length(qg))
    for (a in seq_along(qg)) {
      loc <- nn$nn.idx[a, ]; dd <- nn$nn.dists[a, ]
      global <- data_idx[loc]
      keep <- global != qg[a] & is.finite(dd)
      global <- global[keep]; dd <- dd[keep]
      if (length(global) > k_use) {
        global <- global[seq_len(k_use)]; dd <- dd[seq_len(k_use)]
      }
      cand[[a]] <- global; dist[[a]] <- dd
      if (length(dd)) {
        posd <- dd[dd > 0]
        h <- if (length(posd)) posd[min(length(posd), ctl$local_k)] else 1
        if (!is.finite(h) || h <= 0) h <- 1
        hvec[a] <- h
      } else hvec[a] <- 1
    }
    pair_n <- vapply(cand, length, integer(1))
    if (!sum(pair_n)) next
    q_flat <- rep(qg, pair_n)
    d_flat <- unlist(cand, use.names = FALSE)
    de_flat <- unlist(dist, use.names = FALSE)
    h_flat <- rep(hvec, pair_n)
    qm <- membership[q_flat]; dm <- membership[d_flat]
    tree_index <- .dk_tree_index_for_local_stratum(membership_fit, s)
    dt_flat <- .dk_membership_tree_pair_distance(tree_index, cells[data_idx], membership[data_idx], qm, dm)
    tau <- tree_tau
    finite_tree <- dt_flat[is.finite(dt_flat) & dt_flat > 0]
    if (is.null(tau)) tau <- if (length(finite_tree)) stats::median(finite_tree) else 1
    if (!is.finite(tau) || tau <= 0) tau <- 1
    offset <- cumsum(c(0L, pair_n))
    for (a in seq_along(qg)) {
      if (!pair_n[a]) next
      idx <- seq.int(offset[a] + 1L, offset[a + 1L])
      donors <- d_flat[idx]; de <- de_flat[idx]; dt <- dt_flat[idx]; h <- hvec[a]
      emb_term <- de^2 / (2 * h^2)
      cost <- emb_term
      have_tree <- is.finite(dt)
      cost[have_tree] <- tree_weight * (dt[have_tree] / tau) + (1 - tree_weight) * emb_term[have_tree]
      cost <- cost + ctl$membership_penalty * as.numeric(membership[donors] != membership[qg[a]])
      logp <- -cost
      logp <- logp - max(logp)
      prior <- exp(logp)
      if (!all(is.finite(prior)) || sum(prior) <= 0) prior <- rep(1, length(donors))
      prior <- prior / sum(prior)
      delta <- sweep(z[donors, , drop = FALSE], 2L, z[qg[a], ], "-") / h
      fit <- .dk_barycentric_weights(delta, prior, ctl$barycentric_lambda,
                                     ctl$barycentric_iter, ctl$barycentric_tol)
      w <- fit$weight
      qr <- unname(qrow[as.character(qg[a])])
      bandwidth[qr] <- h; candidate_count[qr] <- length(donors)
      state_error[qr] <- fit$state_error; prior_error[qr] <- fit$prior_error
      iterations[qr] <- fit$iterations; converged[qr] <- fit$converged
      embed_mean[qr] <- sum(w * de)
      tree_ok <- is.finite(dt)
      if (any(tree_ok) && sum(w[tree_ok]) > 0)
        tree_mean[qr] <- sum(w[tree_ok] * dt[tree_ok]) / sum(w[tree_ok])
      same_membership_weight[qr] <- sum(w[membership[donors] == membership[qg[a]]])
      wi <- c(wi, rep.int(qr, length(donors)))
      wj <- c(wj, donors)
      wx <- c(wx, w)
      px <- c(px, prior)
    }
  }
  W <- Matrix::sparseMatrix(i = wi, j = wj, x = wx, dims = c(nq, n),
                            dimnames = list(cells[query_cells], cells))
  P <- Matrix::sparseMatrix(i = wi, j = wj, x = px, dims = c(nq, n),
                            dimnames = list(cells[query_cells], cells))
  list(W = W, W2 = W * W, prior = P, prior2 = P * P, query_cells = query_cells,
       bandwidth = bandwidth, candidate_count = candidate_count,
       state_error = state_error, prior_error = prior_error,
       tree_distance_weighted_mean = tree_mean,
       embedding_distance_weighted_mean = embed_mean,
       same_membership_weight = same_membership_weight,
       iterations = iterations, converged = converged,
       cell_stratum = as.character(strata), control = ctl)
}

.dk_barycentric_stats_batch <- function(x, geometry) {
  W <- geometry$W; W2 <- geometry$W2
  nq <- nrow(W); ng <- nrow(x)
  if (!nq || !ng) {
    empty <- matrix(numeric(), ng, nq)
    return(list(mean = empty, variance = empty, effective_n = empty,
                prevalence = empty, n_donors = empty))
  }
  Wt <- Matrix::t(W); W2t <- Matrix::t(W2)
  pos <- x > 0
  num <- as.matrix(x %*% Wt)
  den <- as.matrix(pos %*% Wt)
  ss <- as.matrix((x * x) %*% Wt)
  den2 <- as.matrix(pos %*% W2t)
  mu <- matrix(NA_real_, ng, nq); ok <- den > 0; mu[ok] <- num[ok] / den[ok]
  neff <- matrix(0, ng, nq); nok <- den2 > 0; neff[nok] <- den[nok]^2 / den2[nok]
  var_num <- ss - 2 * mu * num + mu * mu * den
  var_den <- den - den2 / pmax(den, .Machine$double.eps)
  var <- matrix(NA_real_, ng, nq); vok <- is.finite(var_num) & is.finite(var_den) & var_den > 0
  var[vok] <- pmax(var_num[vok] / var_den[vok], 0)
  total <- as.numeric(Matrix::rowSums(W))
  prevalence <- matrix(NA_real_, ng, nq); tok <- total > 0
  if (any(tok)) prevalence[, tok] <- sweep(den[, tok, drop = FALSE], 2L, total[tok], "/")
  support <- W
  support@x[] <- 1
  n_donors <- as.matrix((pos * 1) %*% Matrix::t(support))
  list(mean = mu, variance = var, effective_n = neff,
       prevalence = prevalence, n_donors = n_donors)
}

.dk_barycentric_predict_events <- function(x, embedding, membership, events,
                                           membership_fit = NULL, hard_stratum = NULL,
                                           tree_weight = 0.5, tree_tau = NULL,
                                           local_k = 30L, candidate_k = 100L,
                                           membership_penalty = 1,
                                           barycentric_lambda = 1,
                                           barycentric_iter = 20L,
                                           barycentric_tol = 1e-6,
                                           min_effective_donors = 5) {
  n_ev <- nrow(events)
  out <- list(
    prediction = rep(NA_real_, n_ev), factor_prediction = rep(NA_real_, n_ev),
    prediction_sd = rep(NA_real_, n_ev), predictability = rep(NA_real_, n_ev),
    shrinkage = rep(NA_real_, n_ev), n_observed_gene = integer(n_ev),
    factor_rank = integer(n_ev), factor_features = integer(n_ev),
    factor_iterations = integer(n_ev), factor_converged = logical(n_ev),
    recovery_method = rep("unavailable", n_ev), target_mode = rep("positive", n_ev),
    cell_prediction = rep(NA_real_, n_ev), cell_available = rep(FALSE, n_ev),
    n_donors = integer(n_ev), bandwidth = rep(NA_real_, n_ev),
    local_positive_mean = rep(NA_real_, n_ev), local_positive_variance = rep(NA_real_, n_ev),
    local_positive_prevalence = rep(NA_real_, n_ev), effective_donors = numeric(n_ev),
    tree_distance_weighted_mean = rep(NA_real_, n_ev),
    embedding_distance_weighted_mean = rep(NA_real_, n_ev)
  )
  if (!n_ev) return(c(out, list(geometry = NULL)))
  geometry <- .dk_build_barycentric_geometry(
    embedding, membership, unique(events$j), membership_fit, hard_stratum,
    tree_weight, tree_tau, local_k, candidate_k, membership_penalty,
    barycentric_lambda, barycentric_iter, barycentric_tol, min_effective_donors
  )
  qpos <- match(events$j, geometry$query_cells)
  gene_index <- split(seq_len(n_ev), events$i)
  genes <- as.integer(names(gene_index))
  nq <- length(geometry$query_cells)
  batch_size <- max(1L, min(128L, as.integer(floor(2e6 / max(1L, nq)))))
  for (start in seq.int(1L, length(genes), by = batch_size)) {
    gset <- genes[start:min(length(genes), start + batch_size - 1L)]
    sb <- .dk_barycentric_stats_batch(x[gset, , drop = FALSE], geometry)
    for (kk in seq_along(gset)) {
      g <- gset[kk]; qg <- gene_index[[as.character(g)]]; qp <- qpos[qg]
      mu <- sb$mean[kk, qp]; vv <- sb$variance[kk, qp]
      ne <- sb$effective_n[kk, qp]; nd <- as.integer(round(sb$n_donors[kk, qp]))
      supported <- is.finite(mu) & is.finite(ne) & ne >= min_effective_donors
      pred <- rep(NA_real_, length(qg)); pred[supported] <- pmax(mu[supported], 0)
      sd <- rep(NA_real_, length(qg))
      vok <- supported & is.finite(vv) & vv >= 0
      sd[vok] <- sqrt(vv[vok] * (1 + 1 / pmax(ne[vok], 1)))
      pe <- geometry$prior_error[qp]; se <- geometry$state_error[qp]
      gain <- rep(NA_real_, length(qg)); gok <- is.finite(pe) & pe > 0 & is.finite(se)
      gain[gok] <- pmax(0, pmin(1, 1 - se[gok] / pe[gok]))
      out$prediction[qg] <- pred; out$prediction_sd[qg] <- sd
      out$predictability[qg] <- gain; out$n_observed_gene[qg] <- nd
      out$recovery_method[qg[supported]] <- "barycentric_mean"
      out$cell_prediction[qg] <- pred; out$cell_available[qg] <- supported
      out$n_donors[qg] <- nd; out$bandwidth[qg] <- geometry$bandwidth[qp]
      out$local_positive_mean[qg] <- mu; out$local_positive_variance[qg] <- vv
      out$local_positive_prevalence[qg] <- sb$prevalence[kk, qp]
      out$effective_donors[qg] <- ne
      out$tree_distance_weighted_mean[qg] <- geometry$tree_distance_weighted_mean[qp]
      out$embedding_distance_weighted_mean[qg] <- geometry$embedding_distance_weighted_mean[qp]
    }
  }
  c(out, list(geometry = geometry))
}
