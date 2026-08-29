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
  if (!is.null(hard_stratum)) {
    s <- .dk_align_vector(hard_stratum, cells, "hard_stratum")
    s <- as.character(s)
    return(factor(s, levels = unique(s)))
  }
  has_hard <- !is.null(membership_fit) && isTRUE(membership_fit$settings$has_hard_stratum)
  if (has_hard && !is.null(membership_fit$cell_stratum)) {
    s <- membership_fit$cell_stratum
    if (!is.null(names(s))) s <- s[cells]
    if (length(s) == length(cells) && !anyNA(s)) return(factor(as.character(s), levels = unique(as.character(s))))
  }
  factor(paste0("membership_", membership), levels = unique(paste0("membership_", membership)))
}

.dk_tree_index_for_local_stratum <- function(membership_fit, stratum) {
  if (is.null(membership_fit) || !length(membership_fit$tree_indices)) return(NULL)
  ti <- membership_fit$tree_indices
  if (!is.null(names(ti)) && stratum %in% names(ti)) return(ti[[stratum]])
  if (length(ti) == 1L) return(ti[[1L]])
  NULL
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
  block <- interaction(strata, factor(membership, levels = unique(membership)), drop = TRUE, lex.order = TRUE)
  block_id <- as.integer(block)
  i_chunks <- list(); j_chunks <- list(); w_chunks <- list(); dt_chunks <- list(); de_chunks <- list(); nchunk <- 0L
  bandwidth <- rep(NA_real_, n); candidate_count <- integer(n)
  for (bid in seq_len(nlevels(block))) {
    idx <- which(block_id == bid); nb <- length(idx)
    if (nb <= 1L) next
    zb <- z[idx, , drop = FALSE]
    k_use <- min(nb - 1L, ctl$candidate_k)
    nn <- RANN::nn2(data = zb, query = zb, k = min(nb, k_use + 1L))
    s <- as.character(strata[idx[1L]])
    tree_index <- .dk_tree_index_for_local_stratum(membership_fit, s)
    cand_list <- vector("list", nb); de_list <- vector("list", nb); dt_list <- vector("list", nb)
    for (a in seq_len(nb)) {
      cand <- nn$nn.idx[a, ]; cand <- cand[cand != a]
      if (length(cand) > k_use) cand <- cand[seq_len(k_use)]
      if (!length(cand)) next
      delta <- zb[cand, , drop = FALSE] - matrix(zb[a, ], nrow = length(cand), ncol = ncol(zb), byrow = TRUE)
      de <- sqrt(rowSums(delta * delta))
      ord <- order(de, cand); cand <- cand[ord]; de <- de[ord]
      pos_de <- de[is.finite(de) & de > 0]
      h <- if (length(pos_de)) pos_de[min(length(pos_de), ctl$local_k)] else 1
      if (!is.finite(h) || h <= 0) h <- 1
      bandwidth[idx[a]] <- h; candidate_count[idx[a]] <- length(cand)
      dt <- .dk_tree_distance(tree_index, rep(cells[idx[a]], length(cand)), cells[idx[cand]])
      cand_list[[a]] <- cand; de_list[[a]] <- de; dt_list[[a]] <- dt
    }
    finite_tree <- unlist(lapply(dt_list, function(v) v[is.finite(v) & v > 0]), use.names = FALSE)
    tau <- ctl$tree_tau
    if (is.null(tau)) tau <- if (length(finite_tree)) stats::median(finite_tree) else 1
    if (!is.finite(tau) || tau <= 0) tau <- 1
    for (a in seq_len(nb)) {
      cand <- cand_list[[a]]; if (!length(cand)) next
      de <- de_list[[a]]; dt <- dt_list[[a]]; h <- bandwidth[idx[a]]
      have_tree <- is.finite(dt)
      alpha <- ifelse(have_tree, ctl$tree_weight, 0)
      tree_term <- numeric(length(cand)); tree_term[have_tree] <- dt[have_tree] / tau
      emb_term <- de^2 / (2 * h^2)
      w <- exp(-alpha * tree_term - (1 - alpha) * emb_term)
      keep <- is.finite(w) & w > 1e-12
      if (!any(keep)) next
      nchunk <- nchunk + 1L
      i_chunks[[nchunk]] <- rep.int(idx[a], sum(keep)); j_chunks[[nchunk]] <- idx[cand[keep]]
      w_chunks[[nchunk]] <- w[keep]
      dt_keep <- dt[keep]; dt_keep[!is.finite(dt_keep)] <- NA_real_
      dt_chunks[[nchunk]] <- dt_keep; de_chunks[[nchunk]] <- de[keep]
    }
  }
  if (nchunk) {
    ii <- unlist(i_chunks, use.names = FALSE); jj <- unlist(j_chunks, use.names = FALSE)
    ww <- unlist(w_chunks, use.names = FALSE); dt_all <- unlist(dt_chunks, use.names = FALSE)
    de_all <- unlist(de_chunks, use.names = FALSE)
  } else {
    ii <- integer(); jj <- integer(); ww <- numeric(); dt_all <- numeric(); de_all <- numeric()
  }
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n), dimnames = list(cells, cells))
  W2 <- W * W
  total_weight <- as.numeric(Matrix::rowSums(W))
  tree_ok <- is.finite(dt_all)
  Dtree <- Matrix::sparseMatrix(i = ii[tree_ok], j = jj[tree_ok], x = dt_all[tree_ok],
                                dims = c(n, n), dimnames = list(cells, cells))
  Dembed <- Matrix::sparseMatrix(i = ii, j = jj, x = de_all, dims = c(n, n), dimnames = list(cells, cells))
  embed_mean <- rep(NA_real_, n); ok <- total_weight > 0
  if (any(ok)) embed_mean[ok] <- as.numeric(Matrix::rowSums(W * Dembed))[ok] / total_weight[ok]
  tree_mean <- rep(NA_real_, n)
  if (length(tree_ok) && any(tree_ok)) {
    tree_mask <- Dtree != 0; tree_w <- W * tree_mask
    tree_den <- as.numeric(Matrix::rowSums(tree_w)); tok <- tree_den > 0
    if (any(tok)) tree_mean[tok] <- as.numeric(Matrix::rowSums(W * Dtree))[tok] / tree_den[tok]
  }
  list(W = W, W2 = W2, total_weight = total_weight,
       tree_distance = Dtree, embedding_distance = Dembed,
       tree_distance_weighted_mean = tree_mean,
       embedding_distance_weighted_mean = embed_mean,
       bandwidth = bandwidth, candidate_count = candidate_count,
       cell_stratum = as.character(strata), block_id = block_id, control = ctl)
}

.dk_local_gene_stats <- function(xg, geometry) {
  W <- geometry$W
  W2 <- geometry$W2; if (is.null(W2)) W2 <- W * W
  total_w <- geometry$total_weight; if (is.null(total_w)) total_w <- as.numeric(Matrix::rowSums(W))
  pos <- is.finite(xg) & xg > 0
  xp <- numeric(length(xg)); xp[pos] <- xg[pos]
  den <- as.numeric(W %*% as.numeric(pos)); num <- as.numeric(W %*% xp)
  ss <- as.numeric(W %*% (xp * xp)); den2 <- as.numeric(W2 %*% as.numeric(pos))
  mu <- rep(NA_real_, length(xg)); ok <- den > 0; mu[ok] <- num[ok] / den[ok]
  neff <- rep(0, length(xg)); neff[den2 > 0] <- den[den2 > 0]^2 / den2[den2 > 0]
  var <- rep(NA_real_, length(xg))
  var_num <- ss - 2 * mu * num + mu * mu * den
  var_den <- den - den2 / pmax(den, .Machine$double.eps)
  vok <- is.finite(var_num) & is.finite(var_den) & var_den > 0
  var[vok] <- pmax(var_num[vok] / var_den[vok], 0)
  prevalence <- rep(NA_real_, length(xg)); pok <- total_w > 0; prevalence[pok] <- den[pok] / total_w[pok]
  list(mean = mu, variance = var, effective_n = neff, prevalence = prevalence,
       positive_weight = den, total_weight = total_w, positive = pos)
}

.dk_local_gene_stats_batch <- function(x, geometry) {
  W <- geometry$W; W2 <- geometry$W2; if (is.null(W2)) W2 <- W * W
  total_w <- geometry$total_weight; if (is.null(total_w)) total_w <- as.numeric(Matrix::rowSums(W))
  pos <- x > 0
  Wt <- geometry$Wt; if (is.null(Wt)) Wt <- Matrix::t(W)
  W2t <- geometry$W2t; if (is.null(W2t)) W2t <- Matrix::t(W2)
  num <- as.matrix(x %*% Wt); den <- as.matrix(pos %*% Wt)
  ss <- as.matrix((x * x) %*% Wt); den2 <- as.matrix(pos %*% W2t)
  mu <- matrix(NA_real_, nrow(x), ncol(x)); ok <- den > 0; mu[ok] <- num[ok] / den[ok]
  neff <- matrix(0, nrow(x), ncol(x)); nok <- den2 > 0; neff[nok] <- den[nok]^2 / den2[nok]
  var_num <- ss - 2 * mu * num + mu * mu * den
  var_den <- den - den2 / pmax(den, .Machine$double.eps)
  var <- matrix(NA_real_, nrow(x), ncol(x)); vok <- is.finite(var_num) & is.finite(var_den) & var_den > 0
  var[vok] <- pmax(var_num[vok] / var_den[vok], 0)
  prevalence <- matrix(NA_real_, nrow(x), ncol(x)); pok <- total_w > 0
  if (any(pok)) prevalence[, pok] <- sweep(den[, pok, drop = FALSE], 2L, total_w[pok], "/")
  list(mean = mu, variance = var, effective_n = neff, prevalence = prevalence, positive = as.matrix(pos))
}

.dk_weighted_local_residual_target <- function(xg, stats_g, scores, W, query,
                                               ridge = 1, min_target_observed = 20L,
                                               min_effective_donors = 5,
                                               local_info_kappa = 5,
                                               support_adaptive_rank = FALSE) {
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0)
    stop("factor_ridge must be >= 0", call. = FALSE)
  nq <- length(query)
  ans <- list(prediction = rep(NA_real_, nq), residual_prediction = rep(NA_real_, nq),
              prediction_sd = rep(NA_real_, nq), predictability = numeric(nq),
              shrinkage = numeric(nq), effective_n = numeric(nq),
              n_donors = integer(nq), method = rep("unavailable", nq))
  if (!nq) return(ans)
  mu_q <- stats_g$mean[query]; var_q <- stats_g$variance[query]; neff_q <- stats_g$effective_n[query]
  supported <- is.finite(mu_q) & is.finite(neff_q) & neff_q >= min_effective_donors
  base_var <- rep(0, nq); bok <- is.finite(var_q) & neff_q > 0; base_var[bok] <- var_q[bok] / neff_q[bok]
  fallback_var <- base_var; vok <- is.finite(var_q); fallback_var[vok] <- pmax(0, var_q[vok] + base_var[vok])
  if (any(supported)) {
    ans$prediction[supported] <- pmax(0, mu_q[supported])
    ans$prediction_sd[supported] <- sqrt(pmax(0, fallback_var[supported]))
    ans$method[supported] <- "tree_local_mean"
  }
  donor_all <- which(stats_g$positive & is.finite(stats_g$mean))
  if (!any(supported) || !length(donor_all) || is.null(scores) || !ncol(scores)) return(ans)
  residual <- xg - stats_g$mean
  qfit <- query[supported]
  Wq0 <- as.matrix(W[qfit, donor_all, drop = FALSE])
  aggregate_w <- colSums(Wq0)
  keep <- is.finite(aggregate_w) & aggregate_w > 0 & is.finite(residual[donor_all])
  donor <- donor_all[keep]; aggregate_w <- aggregate_w[keep]
  if (length(donor) < max(as.integer(min_target_observed), 2L)) return(ans)
  scores_use <- scores
  if (support_adaptive_rank) {
    k_use <- min(ncol(scores), max(0L, floor((length(donor) - 3L) / 2L)))
    if (k_use < 1L) return(ans)
    scores_use <- scores[, seq_len(k_use), drop = FALSE]
  }
  Xo <- cbind(1, scores_use[donor, , drop = FALSE]); y <- residual[donor]; p <- ncol(Xo)
  if (length(donor) < p + 2L) return(ans)
  scale_w <- aggregate_w / mean(aggregate_w); sw <- sqrt(scale_w)
  Xw <- Xo * sw; yw <- y * sw
  P <- diag(c(0, rep(ridge, p - 1L)), nrow = p)
  A <- crossprod(Xw) + P
  inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(inv) || any(!is.finite(inv))) return(ans)
  beta <- as.vector(inv %*% crossprod(Xw, yw)); fitted <- as.vector(Xo %*% beta)
  hdiag <- scale_w * rowSums((Xo %*% inv) * Xo)
  loo <- y - (y - fitted) / pmax(1 - hdiag, 1e-6)
  den <- sum(aggregate_w * loo * loo)
  qpred <- if (is.finite(den) && den > 1e-12) sum(aggregate_w * loo * y) / den else 0
  qpred <- max(0, min(1, qpred))
  qinfo <- if (local_info_kappa > 0) neff_q / (neff_q + local_info_kappa) else rep(1, nq)
  qfinal <- qpred * qinfo; qfinal[!supported] <- 0
  Xq <- cbind(1, scores_use[query, , drop = FALSE]); raw <- as.vector(Xq %*% beta)
  pred <- pmax(mu_q + qfinal * raw, 0)
  Wqd <- as.matrix(W[query, donor, drop = FALSE])
  ans$n_donors <- rowSums(Wqd > 0); ans$effective_n <- neff_q
  null_sse <- as.numeric(Wqd %*% (y * y)); cross_term <- as.numeric(Wqd %*% (y * loo))
  loo_sse <- as.numeric(Wqd %*% (loo * loo))
  final_sse <- pmax(null_sse - 2 * qfinal * cross_term + qfinal^2 * loo_sse, 0)
  wsum <- rowSums(Wqd)
  predability <- numeric(nq); pok <- supported & is.finite(null_sse) & null_sse > 0
  predability[pok] <- pmax(0, pmin(1, 1 - final_sse[pok] / null_sse[pok]))
  sigma2 <- rep(NA_real_, nq); sok <- supported & wsum > 0; sigma2[sok] <- final_sse[sok] / wsum[sok]
  bad_sigma <- !is.finite(sigma2) | sigma2 < 0; sigma2[bad_sigma] <- var_q[bad_sigma]
  meat <- crossprod(Xo, Xo * scale_w^2); M <- inv %*% meat %*% inv
  hq <- rowSums((Xq %*% M) * Xq)
  pv <- sigma2 * (1 + qfinal^2 * pmax(hq, 0)) + base_var
  factor_use <- supported & qfinal > 0 & is.finite(pred)
  ans$prediction[supported] <- pred[supported]
  ans$residual_prediction[supported] <- raw[supported]
  ans$predictability[supported] <- predability[supported]
  ans$shrinkage[supported] <- qfinal[supported]
  if (any(factor_use)) {
    ans$prediction_sd[factor_use] <- sqrt(pmax(0, pv[factor_use]))
    ans$method[factor_use] <- "tree_local_factor"
  }
  ans$factor_rank_used <- ncol(scores_use)
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
                                          local_info_kappa = 5,
                                          factor_crossfit_folds = 1L,
                                          factor_crossfit_seed = 1L,
                                          support_adaptive_rank = FALSE) {
  factor_crossfit_folds <- as.integer(factor_crossfit_folds)[1L]
  factor_crossfit_seed <- as.integer(factor_crossfit_seed)[1L]
  if (!is.finite(factor_crossfit_folds) || factor_crossfit_folds < 1L)
    stop("factor_crossfit_folds must be >= 1", call. = FALSE)
  if (!is.finite(factor_crossfit_seed))
    stop("factor_crossfit_seed must be finite", call. = FALSE)
  n_ev <- nrow(events)
  out <- list(
    prediction = rep(NA_real_, n_ev), factor_prediction = rep(NA_real_, n_ev),
    prediction_sd = rep(NA_real_, n_ev), predictability = numeric(n_ev),
    shrinkage = numeric(n_ev), n_observed_gene = integer(n_ev),
    factor_rank = integer(n_ev), factor_features = integer(n_ev),
    factor_fold = integer(n_ev), bias_calibration = numeric(n_ev),
    factor_iterations = integer(n_ev), factor_converged = logical(n_ev),
    recovery_method = rep("unavailable", n_ev), target_mode = rep("positive", n_ev),
    local_positive_mean = rep(NA_real_, n_ev), local_positive_variance = rep(NA_real_, n_ev),
    local_positive_prevalence = rep(NA_real_, n_ev), effective_donors = numeric(n_ev),
    tree_distance_weighted_mean = rep(NA_real_, n_ev),
    embedding_distance_weighted_mean = rep(NA_real_, n_ev)
  )
  if (!n_ev) return(c(out, list(geometry = NULL)))
  geometry <- .dk_build_local_geometry(
    embedding, membership, membership_fit, hard_stratum,
    tree_weight, tree_tau, local_k, candidate_k,
    min_effective_donors, local_info_kappa
  )
  event_block <- geometry$block_id[events$j]
  event_by_block <- split(seq_len(n_ev), event_block)
  batch_size <- 256L
  for (b in names(event_by_block)) {
    bid <- as.integer(b); qev <- event_by_block[[b]]; cells <- which(geometry$block_id == bid)
    ev_i <- events$i[qev]; ev_j <- events$j[qev]
    local_events <- data.frame(i = ev_i, j = match(ev_j, cells))
    Wb <- geometry$W[cells, cells, drop = FALSE]; W2b <- geometry$W2[cells, cells, drop = FALSE]
    Wb_dense <- as.matrix(Wb)
    geob <- list(W = Wb, W2 = W2b, Wt = Matrix::t(Wb), W2t = Matrix::t(W2b), total_weight = geometry$total_weight[cells])
    gene_index <- split(seq_along(ev_i), ev_i)
    genes <- as.integer(names(gene_index))
    folds <- .dk_target_fold(genes, factor_crossfit_folds, factor_crossfit_seed)
    for (fold in unique(folds)) {
      fold_genes <- genes[folds == fold]
      fold_event <- ev_i %in% fold_genes
      fit <- .dk_membership_factor_scores(
        x, cells, local_events[fold_event, , drop = FALSE],
        rank = factor_rank, feature_max = factor_features,
        min_feature_observed = min_feature_observed
      )
      scores <- if (is.null(fit)) NULL else fit$scores
      for (start in seq.int(1L, length(fold_genes), by = batch_size)) {
        gset <- fold_genes[start:min(length(fold_genes), start + batch_size - 1L)]
        xb <- x[gset, cells, drop = FALSE]; sb <- .dk_local_gene_stats_batch(xb, geob)
        qg0_batch <- which(ev_i %in% gset)
        qg_batch <- qev[qg0_batch]
        query_gene <- match(ev_i[qg0_batch], gset)
        query_global <- ev_j[qg0_batch]
        query <- match(query_global, cells)
        scores_cpp <- if (is.null(scores)) matrix(numeric(), length(cells), 0L) else scores
        cpp <- dk_batch_tree_ridge_positive_cpp(
          as.matrix(xb), sb$mean, sb$variance, sb$effective_n,
          Wb_dense, scores_cpp, query_gene, query,
          factor_ridge, as.integer(min_target_observed),
          min_effective_donors, local_info_kappa,
          support_adaptive_rank
        )
        out$prediction[qg_batch] <- cpp$prediction
        out$factor_prediction[qg_batch] <- cpp$residual_prediction
        out$prediction_sd[qg_batch] <- cpp$prediction_sd
        out$predictability[qg_batch] <- cpp$predictability
        out$shrinkage[qg_batch] <- cpp$shrinkage
        out$n_observed_gene[qg_batch] <- cpp$n_donors
        out$recovery_method[qg_batch] <- c(
          "unavailable", "tree_local_mean", "tree_local_factor"
        )[cpp$method + 1L]
        out$factor_rank[qg_batch] <- cpp$factor_rank_used
        out$factor_fold[qg_batch] <- fold
        stat_idx <- cbind(query_gene, query)
        out$local_positive_mean[qg_batch] <- sb$mean[stat_idx]
        out$local_positive_variance[qg_batch] <- sb$variance[stat_idx]
        out$local_positive_prevalence[qg_batch] <- sb$prevalence[stat_idx]
        out$effective_donors[qg_batch] <- sb$effective_n[stat_idx]
        out$tree_distance_weighted_mean[qg_batch] <- geometry$tree_distance_weighted_mean[query_global]
        out$embedding_distance_weighted_mean[qg_batch] <- geometry$embedding_distance_weighted_mean[query_global]
        if (!is.null(fit)) {
          out$factor_features[qg_batch] <- fit$n_features
          out$factor_iterations[qg_batch] <- fit$iterations
          out$factor_converged[qg_batch] <- fit$converged
        }
      }
    }
  }
  c(out, list(geometry = geometry))
}
