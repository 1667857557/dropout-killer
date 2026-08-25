.dk_tree_local_control <- function(tree_weight = 0.5, tree_tau = NULL, local_k = 30L,
                                   candidate_k = 100L, min_effective_donors = 5,
                                   local_info_kappa = 5) {
  if (!is.numeric(tree_weight) || length(tree_weight) != 1L || !is.finite(tree_weight) || tree_weight < 0 || tree_weight > 1)
    stop("tree_weight must be in [0,1]", call. = FALSE)
  if (!is.null(tree_tau) && (!is.numeric(tree_tau) || length(tree_tau) != 1L || !is.finite(tree_tau) || tree_tau <= 0))
    stop("tree_tau must be NULL or a positive finite number", call. = FALSE)
  for (z in list(local_k = local_k, candidate_k = candidate_k, min_effective_donors = min_effective_donors)) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z < 1) stop("local neighborhood controls must be >= 1", call. = FALSE)
  }
  if (!is.numeric(local_info_kappa) || length(local_info_kappa) != 1L || !is.finite(local_info_kappa) || local_info_kappa < 0)
    stop("local_info_kappa must be >= 0", call. = FALSE)
  list(tree_weight = as.numeric(tree_weight), tree_tau = tree_tau,
       local_k = as.integer(local_k), candidate_k = as.integer(candidate_k),
       min_effective_donors = as.numeric(min_effective_donors),
       local_info_kappa = as.numeric(local_info_kappa))
}

.dk_resolve_local_stratum <- function(cells, membership, membership_fit = NULL, hard_stratum = NULL) {
  if (!is.null(hard_stratum)) return(.dk_align_vector(hard_stratum, cells, "hard_stratum"))
  if (!is.null(membership_fit) && !is.null(membership_fit$cell_stratum)) {
    s <- membership_fit$cell_stratum
    if (!is.null(names(s))) s <- s[cells]
    if (length(s) == length(cells) && !anyNA(s)) return(factor(as.character(s), levels = unique(as.character(s))))
  }
  factor(paste0("membership_", membership), levels = unique(paste0("membership_", membership)))
}

.dk_build_local_geometry <- function(embedding, membership, membership_fit = NULL,
                                     hard_stratum = NULL, tree_weight = 0.5,
                                     tree_tau = NULL, local_k = 30L,
                                     candidate_k = 100L,
                                     min_effective_donors = 5,
                                     local_info_kappa = 5) {
  ctl <- .dk_tree_local_control(tree_weight, tree_tau, local_k, candidate_k,
                                min_effective_donors, local_info_kappa)
  z <- as.matrix(embedding); n <- nrow(z); cells <- rownames(z)
  if (is.null(cells)) stop("embedding must have cell row names for tree-local recovery", call. = FALSE)
  strata <- .dk_resolve_local_stratum(cells, membership, membership_fit, hard_stratum)
  ii <- integer(); jj <- integer(); ww <- numeric(); dt_all <- numeric(); de_all <- numeric()
  bandwidth <- rep(NA_real_, n); candidate_count <- integer(n)
  lev <- levels(strata)
  for (s in lev) {
    idx <- which(strata == s); ns <- length(idx)
    if (ns <= 1L) next
    zs <- z[idx, , drop = FALSE]
    kq <- min(ns, ctl$candidate_k + 1L)
    nn <- RANN::nn2(data = zs, query = zs, k = kq)
    tree_index <- NULL
    if (!is.null(membership_fit) && length(membership_fit$tree_indices)) tree_index <- membership_fit$tree_indices[[s]]
    cand_list <- vector("list", ns); de_list <- vector("list", ns); dt_list <- vector("list", ns)
    for (a in seq_len(ns)) {
      kn <- nn$nn.idx[a, ]; kn <- kn[kn != a]
      same_m <- which(membership[idx] == membership[idx[a]])
      cand <- sort(unique(c(kn, same_m)))
      cand <- cand[cand != a]
      if (!length(cand)) next
      de <- sqrt(rowSums((zs[cand, , drop = FALSE] - matrix(zs[a, ], nrow = length(cand), ncol = ncol(zs), byrow = TRUE))^2))
      ord <- order(de, cand)
      cand <- cand[ord]; de <- de[ord]
      pos_de <- de[is.finite(de) & de > 0]
      h <- if (length(pos_de)) pos_de[min(length(pos_de), ctl$local_k)] else 1
      if (!is.finite(h) || h <= 0) h <- 1
      bandwidth[idx[a]] <- h
      candidate_count[idx[a]] <- length(cand)
      qname <- rep(cells[idx[a]], length(cand)); dname <- cells[idx[cand]]
      dt <- .dk_tree_distance(tree_index, qname, dname)
      cand_list[[a]] <- cand; de_list[[a]] <- de; dt_list[[a]] <- dt
    }
    finite_tree <- unlist(lapply(dt_list, function(v) v[is.finite(v) & v > 0]), use.names = FALSE)
    tau <- ctl$tree_tau
    if (is.null(tau)) tau <- if (length(finite_tree)) stats::median(finite_tree) else 1
    if (!is.finite(tau) || tau <= 0) tau <- 1
    for (a in seq_len(ns)) {
      cand <- cand_list[[a]]; if (!length(cand)) next
      de <- de_list[[a]]; dt <- dt_list[[a]]; h <- bandwidth[idx[a]]
      have_tree <- is.finite(dt)
      alpha <- ifelse(have_tree, ctl$tree_weight, 0)
      tree_term <- numeric(length(cand)); tree_term[have_tree] <- dt[have_tree] / tau
      emb_term <- de^2 / (2 * h^2)
      w <- exp(-alpha * tree_term - (1 - alpha) * emb_term)
      keep <- is.finite(w) & w > 1e-12
      if (!any(keep)) next
      q_global <- idx[a]; d_global <- idx[cand[keep]]
      ii <- c(ii, rep.int(q_global, sum(keep))); jj <- c(jj, d_global); ww <- c(ww, w[keep])
      dt_keep <- dt[keep]; dt_keep[!is.finite(dt_keep)] <- NA_real_
      dt_all <- c(dt_all, dt_keep); de_all <- c(de_all, de[keep])
    }
  }
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n), dimnames = list(cells, cells))
  Dtree <- Matrix::sparseMatrix(i = ii[is.finite(dt_all)], j = jj[is.finite(dt_all)], x = dt_all[is.finite(dt_all)],
                                dims = c(n, n), dimnames = list(cells, cells))
  Dembed <- Matrix::sparseMatrix(i = ii, j = jj, x = de_all, dims = c(n, n), dimnames = list(cells, cells))
  list(W = W, tree_distance = Dtree, embedding_distance = Dembed,
       bandwidth = bandwidth, candidate_count = candidate_count,
       cell_stratum = as.character(strata), control = ctl)
}

.dk_local_gene_stats <- function(xg, geometry) {
  W <- geometry$W
  pos <- is.finite(xg) & xg > 0
  xp <- numeric(length(xg)); xp[pos] <- xg[pos]
  den <- as.numeric(W %*% as.numeric(pos))
  num <- as.numeric(W %*% xp)
  ss <- as.numeric(W %*% (xp * xp))
  W2 <- W * W
  den2 <- as.numeric(W2 %*% as.numeric(pos))
  total_w <- as.numeric(Matrix::rowSums(W))
  mu <- rep(NA_real_, length(xg)); ok <- den > 0
  mu[ok] <- num[ok] / den[ok]
  neff <- rep(0, length(xg)); neff[den2 > 0] <- den[den2 > 0]^2 / den2[den2 > 0]
  var <- rep(NA_real_, length(xg))
  var_num <- ss - 2 * mu * num + mu * mu * den
  var_den <- den - den2 / pmax(den, .Machine$double.eps)
  vok <- is.finite(var_num) & is.finite(var_den) & var_den > 0
  var[vok] <- pmax(var_num[vok] / var_den[vok], 0)
  prevalence <- rep(NA_real_, length(xg)); pok <- total_w > 0
  prevalence[pok] <- den[pok] / total_w[pok]
  list(mean = mu, variance = var, effective_n = neff, prevalence = prevalence,
       positive_weight = den, total_weight = total_w, positive = pos)
}

.dk_weighted_local_residual_target <- function(xg, stats_g, scores, W, query,
                                               ridge = 1, min_target_observed = 20L,
                                               min_effective_donors = 5,
                                               local_info_kappa = 5) {
  nq <- length(query)
  ans <- list(prediction = rep(NA_real_, nq), residual_prediction = rep(NA_real_, nq),
              prediction_sd = rep(NA_real_, nq), predictability = numeric(nq),
              shrinkage = numeric(nq), effective_n = numeric(nq),
              n_donors = integer(nq), method = rep("unavailable", nq))
  if (!nq) return(ans)
  donor_all <- which(stats_g$positive & is.finite(stats_g$mean))
  if (!length(donor_all)) return(ans)
  residual <- xg - stats_g$mean
  for (u in seq_along(query)) {
    q <- query[u]
    mu_q <- stats_g$mean[q]; var_q <- stats_g$variance[q]; neff_q <- stats_g$effective_n[q]
    if (!is.finite(mu_q) || neff_q < min_effective_donors) next
    w <- as.numeric(W[q, donor_all, drop = TRUE])
    keep <- is.finite(w) & w > 0 & is.finite(residual[donor_all])
    donor <- donor_all[keep]; w <- w[keep]
    ans$effective_n[u] <- if (length(w) && sum(w * w) > 0) sum(w)^2 / sum(w * w) else 0
    ans$n_donors[u] <- length(donor)
    base_var <- if (is.finite(var_q) && neff_q > 0) var_q / neff_q else 0
    if (length(donor) < max(as.integer(min_target_observed), 2L) || is.null(scores) || !ncol(scores) ||
        ans$effective_n[u] < min_effective_donors) {
      ans$prediction[u] <- max(0, mu_q)
      ans$prediction_sd[u] <- sqrt(max(0, if (is.finite(var_q)) var_q + base_var else base_var))
      ans$method[u] <- "tree_local_mean"
      next
    }
    Xo <- cbind(1, scores[donor, , drop = FALSE]); y <- residual[donor]
    p <- ncol(Xo)
    if (length(donor) < p + 2L) {
      ans$prediction[u] <- max(0, mu_q); ans$prediction_sd[u] <- sqrt(max(0, base_var)); ans$method[u] <- "tree_local_mean"; next
    }
    sw <- sqrt(w / mean(w)); Xw <- Xo * sw; yw <- y * sw
    P <- diag(c(0, rep(ridge, p - 1L)), nrow = p)
    A <- crossprod(Xw) + P
    inv <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(inv) || any(!is.finite(inv))) {
      ans$prediction[u] <- max(0, mu_q); ans$prediction_sd[u] <- sqrt(max(0, base_var)); ans$method[u] <- "tree_local_mean"; next
    }
    beta <- as.vector(inv %*% crossprod(Xw, yw))
    fitted <- as.vector(Xo %*% beta)
    hdiag <- w / mean(w) * rowSums((Xo %*% inv) * Xo)
    loo <- y - (y - fitted) / pmax(1 - hdiag, 1e-6)
    den <- sum(w * loo * loo)
    qpred <- if (is.finite(den) && den > 1e-12) sum(w * loo * y) / den else 0
    qpred <- max(0, min(1, qpred))
    null_sse <- sum(w * y * y); model_sse <- sum(w * (y - qpred * loo)^2)
    predability <- if (is.finite(null_sse) && null_sse > 0) max(0, min(1, 1 - model_sse / null_sse)) else 0
    qinfo <- if (local_info_kappa > 0) ans$effective_n[u] / (ans$effective_n[u] + local_info_kappa) else 1
    qfinal <- qpred * qinfo
    xq <- c(1, scores[q, , drop = TRUE]); raw <- sum(xq * beta)
    pred <- max(0, mu_q + qfinal * raw)
    sigma2 <- if (sum(w) > 0) model_sse / sum(w) else NA_real_
    if (!is.finite(sigma2) || sigma2 < 0) sigma2 <- if (is.finite(var_q)) var_q else 0
    meat <- crossprod(Xo, Xo * (w / mean(w))^2)
    hq <- as.numeric(t(xq) %*% inv %*% meat %*% inv %*% xq)
    pv <- sigma2 * (1 + qfinal^2 * max(0, hq)) + base_var
    ans$prediction[u] <- pred; ans$residual_prediction[u] <- raw
    ans$prediction_sd[u] <- sqrt(max(0, pv)); ans$predictability[u] <- predability
    ans$shrinkage[u] <- qfinal; ans$method[u] <- if (qfinal > 0) "tree_local_factor" else "tree_local_mean"
  }
  ans
}

.dk_tree_local_predict_events <- function(x, embedding, membership, events, membership_fit = NULL,
                                          hard_stratum = NULL, factor_rank = 5L,
                                          factor_features = 2000L, factor_ridge = 1,
                                          min_feature_observed = 20L,
                                          min_target_observed = 20L,
                                          tree_weight = 0.5, tree_tau = NULL,
                                          local_k = 30L, candidate_k = 100L,
                                          min_effective_donors = 5,
                                          local_info_kappa = 5) {
  n_ev <- nrow(events); n <- ncol(x)
  out <- list(prediction = rep(NA_real_, n_ev), factor_prediction = rep(NA_real_, n_ev),
              prediction_sd = rep(NA_real_, n_ev), predictability = numeric(n_ev),
              shrinkage = numeric(n_ev), n_observed_gene = integer(n_ev),
              factor_rank = integer(n_ev), factor_features = integer(n_ev),
              factor_iterations = integer(n_ev), factor_converged = logical(n_ev),
              recovery_method = rep("unavailable", n_ev), target_mode = rep("positive", n_ev),
              local_positive_mean = rep(NA_real_, n_ev), local_positive_variance = rep(NA_real_, n_ev),
              local_positive_prevalence = rep(NA_real_, n_ev), effective_donors = numeric(n_ev),
              tree_distance_weighted_mean = rep(NA_real_, n_ev), embedding_distance_weighted_mean = rep(NA_real_, n_ev))
  if (!n_ev) return(c(out, list(geometry = NULL)))
  geometry <- .dk_build_local_geometry(embedding, membership, membership_fit, hard_stratum,
                                       tree_weight, tree_tau, local_k, candidate_k,
                                       min_effective_donors, local_info_kappa)
  stratum <- geometry$cell_stratum
  events$stratum <- stratum[events$j]
  for (s in unique(events$stratum)) {
    cells <- which(stratum == s); qev <- which(events$stratum == s); evs <- events[qev, , drop = FALSE]
    local_events <- evs; local_events$j <- match(evs$j, cells)
    fit <- .dk_membership_factor_scores(x, cells, local_events, rank = factor_rank,
                                        feature_max = factor_features,
                                        min_feature_observed = min_feature_observed)
    scores <- if (is.null(fit)) NULL else fit$scores
    Ws <- geometry$W[cells, cells, drop = FALSE]
    for (g in unique(evs$i)) {
      qg0 <- which(evs$i == g); qg <- qev[qg0]; query_global <- evs$j[qg0]; query <- match(query_global, cells)
      xg <- as.numeric(x[g, cells]); stats_g <- .dk_local_gene_stats(xg, list(W = Ws))
      tg <- .dk_weighted_local_residual_target(xg, stats_g, scores, Ws, query,
                                               ridge = factor_ridge,
                                               min_target_observed = min_target_observed,
                                               min_effective_donors = min_effective_donors,
                                               local_info_kappa = local_info_kappa)
      out$prediction[qg] <- tg$prediction; out$factor_prediction[qg] <- tg$residual_prediction
      out$prediction_sd[qg] <- tg$prediction_sd; out$predictability[qg] <- tg$predictability
      out$shrinkage[qg] <- tg$shrinkage; out$n_observed_gene[qg] <- tg$n_donors
      out$recovery_method[qg] <- tg$method
      out$local_positive_mean[qg] <- stats_g$mean[query]
      out$local_positive_variance[qg] <- stats_g$variance[query]
      out$local_positive_prevalence[qg] <- stats_g$prevalence[query]
      out$effective_donors[qg] <- tg$effective_n
      for (u in seq_along(query_global)) {
        qcell <- query_global[u]
        wrow <- geometry$W[qcell, , drop = TRUE]; nz <- which(as.numeric(wrow) > 0)
        if (length(nz)) {
          wv <- as.numeric(wrow[nz]); sw <- sum(wv)
          if (sw > 0) {
            td <- as.numeric(geometry$tree_distance[qcell, nz, drop = TRUE])
            ed <- as.numeric(geometry$embedding_distance[qcell, nz, drop = TRUE])
            goodt <- is.finite(td) & td > 0
            if (any(goodt)) out$tree_distance_weighted_mean[qg[u]] <- sum(wv[goodt] * td[goodt]) / sum(wv[goodt])
            goode <- is.finite(ed)
            if (any(goode)) out$embedding_distance_weighted_mean[qg[u]] <- sum(wv[goode] * ed[goode]) / sum(wv[goode])
          }
        }
      }
      if (!is.null(fit)) {
        out$factor_rank[qg] <- fit$rank; out$factor_features[qg] <- fit$n_features
        out$factor_iterations[qg] <- fit$iterations; out$factor_converged[qg] <- fit$converged
      }
    }
  }
  c(out, list(geometry = geometry))
}
