.dk_architecture_methods <- function() c(
  "current_p1_crossfit", "p1_local_residual", "p1_tree_loo_blend",
  "hierarchical_partial_pool", "local_simplex", "p1_stabilized_state",
  "fast_p1_lite", "count_eb_p1_local"
)

.dk_build_architecture_geometry <- function(embedding, membership,
                                            membership_fit = NULL,
                                            hard_stratum = NULL,
                                            candidate_k = 30L,
                                            sibling_memberships = 2L,
                                            tree_weight = 0.75) {
  z <- as.matrix(embedding)
  n <- nrow(z); cells <- rownames(z)
  if (is.null(cells)) stop("embedding must have cell row names", call. = FALSE)
  membership <- .dk_align_membership(membership, cells)
  strata <- as.character(.dk_resolve_local_stratum(
    cells, membership, membership_fit, hard_stratum
  ))
  same <- .dk_build_local_geometry(
    z, membership, membership_fit, hard_stratum,
    tree_weight = tree_weight, local_k = max(30L, candidate_k),
    candidate_k = max(100L, candidate_k), min_effective_donors = 1,
    local_info_kappa = 5
  )
  candidate_k <- as.integer(candidate_k)
  sibling_memberships <- as.integer(sibling_memberships)
  candidate_index <- matrix(0L, n, candidate_k)
  candidate_prior <- matrix(0, n, candidate_k)
  sibling_i <- sibling_j <- integer(); sibling_w <- numeric()
  sibling_map <- list()
  for (s in unique(strata)) {
    ids <- which(strata == s); mids <- unique(membership[ids])
    centroids <- t(vapply(mids, function(m)
      colMeans(z[ids[membership[ids] == m], , drop = FALSE]), numeric(ncol(z))))
    rownames(centroids) <- mids
    representative <- vapply(seq_along(mids), function(k) {
      ii <- ids[membership[ids] == mids[k]]
      dd <- rowSums((z[ii, , drop = FALSE] - matrix(
        centroids[k, ], nrow = length(ii), ncol = ncol(z), byrow = TRUE
      ))^2)
      ii[which.min(dd)]
    }, integer(1))
    cd <- as.matrix(stats::dist(centroids)); diag(cd) <- Inf
    tree_index <- .dk_tree_index_for_local_stratum(membership_fit, s)
    td <- matrix(NA_real_, length(mids), length(mids))
    if (!is.null(tree_index) && length(mids) > 1L) {
      for (a in seq_along(mids)) for (b in seq_along(mids)) if (a != b)
        td[a, b] <- .dk_tree_distance(tree_index, cells[representative[a]],
                                      cells[representative[b]])
    }
    finite_td <- td[is.finite(td) & td > 0]
    tau <- if (length(finite_td)) stats::median(finite_td) else 1
    finite_cd <- cd[is.finite(cd) & cd > 0]
    hc <- if (length(finite_cd)) stats::median(finite_cd) else 1
    cost <- (1 - tree_weight) * cd^2 / (2 * hc^2)
    have_tree <- is.finite(td)
    cost[have_tree] <- tree_weight * td[have_tree] / tau +
      (1 - tree_weight) * cd[have_tree]^2 / (2 * hc^2)
    diag(cost) <- Inf
    sib <- lapply(seq_along(mids), function(a) {
      if (length(mids) <= 1L) integer()
      else mids[head(order(cost[a, ], mids), sibling_memberships)]
    })
    names(sib) <- as.character(mids)
    sibling_map[[s]] <- sib
    for (c in ids) {
      m <- membership[c]
      same_pool <- ids[membership[ids] == m & ids != c]
      sib_ids <- sib[[as.character(m)]]
      sibling_pool <- ids[membership[ids] %in% sib_ids]
      order_pool <- function(pool) {
        if (!length(pool)) return(integer())
        d <- sqrt(rowSums((z[pool, , drop = FALSE] - matrix(
          z[c, ], nrow = length(pool), ncol = ncol(z), byrow = TRUE
        ))^2))
        pool[order(d, pool)]
      }
      same_ord <- order_pool(same_pool)
      take_same <- head(same_ord, candidate_k)
      need <- candidate_k - length(take_same)
      sib_ord <- if (need > 0) head(order_pool(sibling_pool), need) else integer()
      cand <- c(take_same, sib_ord)
      if (!length(cand)) next
      de <- sqrt(rowSums((z[cand, , drop = FALSE] - matrix(
        z[c, ], nrow = length(cand), ncol = ncol(z), byrow = TRUE
      ))^2))
      pos_de <- de[is.finite(de) & de > 0]
      h <- if (length(pos_de)) pos_de[min(length(pos_de), candidate_k)] else 1
      if (!is.finite(h) || h <= 0) h <- 1
      dt <- .dk_tree_distance(tree_index, rep(cells[c], length(cand)), cells[cand])
      alpha <- ifelse(is.finite(dt), tree_weight, 0)
      tree_term <- numeric(length(cand)); tree_term[is.finite(dt)] <- dt[is.finite(dt)] / tau
      w <- exp(-alpha * tree_term - (1 - alpha) * de^2 / (2 * h^2))
      w[!is.finite(w) | w < 0] <- 0
      if (sum(w) <= 0) w[] <- 1
      w <- w / sum(w)
      candidate_index[c, seq_along(cand)] <- cand
      candidate_prior[c, seq_along(cand)] <- w
      is_sibling <- membership[cand] != m
      if (any(is_sibling)) {
        sibling_i <- c(sibling_i, rep.int(c, sum(is_sibling)))
        sibling_j <- c(sibling_j, cand[is_sibling])
        sibling_w <- c(sibling_w, w[is_sibling])
      }
    }
  }
  W_sibling <- Matrix::sparseMatrix(
    i = sibling_i, j = sibling_j, x = sibling_w,
    dims = c(n, n), dimnames = list(cells, cells)
  )
  list(
    same = same, candidate_index = candidate_index,
    candidate_prior = candidate_prior, W_sibling = W_sibling,
    W2_sibling = W_sibling * W_sibling,
    sibling_map = sibling_map, hard_stratum = strata,
    membership = membership,
    control = list(candidate_k = candidate_k,
                   sibling_memberships = sibling_memberships,
                   tree_weight = tree_weight)
  )
}

.dk_positive_group_moments <- function(x, group) {
  lev <- unique(group); G <- nrow(x)
  mean <- variance <- matrix(NA_real_, G, length(lev),
                             dimnames = list(rownames(x), lev))
  n <- matrix(0, G, length(lev), dimnames = list(rownames(x), lev))
  for (k in seq_along(lev)) {
    y <- x[, group == lev[k], drop = FALSE]
    nk <- as.numeric(.dk_row_sums(y > 0)); s1 <- as.numeric(.dk_row_sums(y))
    s2 <- as.numeric(.dk_row_sums(y * y)); ok <- nk > 0
    mean[ok, k] <- s1[ok] / nk[ok]
    vok <- nk > 1
    variance[vok, k] <- pmax((s2[vok] - s1[vok]^2 / nk[vok]) / (nk[vok] - 1), 0)
    n[, k] <- nk
  }
  list(mean = mean, variance = variance, n = n, levels = lev)
}

.dk_local_batch_moments <- function(y, W, W2 = W * W) {
  pos <- y > 0; Wt <- Matrix::t(W); W2t <- Matrix::t(W2)
  den <- as.matrix(pos %*% Wt); den2 <- as.matrix(pos %*% W2t)
  num <- as.matrix(y %*% Wt); ss <- as.matrix((y * y) %*% Wt)
  mu <- matrix(NA_real_, nrow(y), ncol(y)); ok <- den > 0; mu[ok] <- num[ok] / den[ok]
  neff <- matrix(0, nrow(y), ncol(y)); nok <- den2 > 0; neff[nok] <- den[nok]^2 / den2[nok]
  var_num <- ss - 2 * mu * num + mu^2 * den
  var_den <- den - den2 / pmax(den, .Machine$double.eps)
  variance <- matrix(NA_real_, nrow(y), ncol(y)); vok <- is.finite(var_den) & var_den > 0
  variance[vok] <- pmax(var_num[vok] / var_den[vok], 0)
  list(mean = mu, variance = variance, effective_n = neff,
       positive_weight = den, positive_weight2 = den2)
}

.dk_stratum_phi_fallback <- function(counts, size_factor, stratum,
                                     phi_floor = 1e-4) {
  lev <- unique(stratum); out <- matrix(phi_floor, nrow(counts), length(lev),
                                        dimnames = list(rownames(counts), lev))
  for (k in seq_along(lev)) {
    cells <- which(stratum == lev[k]); y <- counts[, cells, drop = FALSE]
    sf <- size_factor[cells]
    u <- y %*% Matrix::Diagonal(x = 1 / sf)
    n <- as.numeric(.dk_row_sums(y > 0)); s1 <- as.numeric(.dk_row_sums(u))
    mu <- s1 / pmax(n, 1)
    centered2 <- as.numeric(.dk_row_sums(u * u)) - 2 * mu * s1 + n * mu^2
    poisson <- as.numeric(.dk_row_sums((y > 0) %*% Matrix::Diagonal(x = 1 / sf))) * mu
    numerator <- pmax(centered2 - poisson, 0)
    phi <- numerator / pmax(n * mu^2, .Machine$double.eps)
    phi[!is.finite(phi)] <- phi_floor
    out[, k] <- pmax(phi, phi_floor)
  }
  out
}

.dk_gamma_power_posterior <- function(prior_mean, phi, weighted_count = 0,
                                      weighted_exposure = 0,
                                      mu_floor = 1e-8,
                                      phi_floor = 1e-4) {
  prior_mean <- pmax(as.numeric(prior_mean), mu_floor)
  phi <- pmax(as.numeric(phi), phi_floor)
  weighted_count <- pmax(as.numeric(weighted_count), 0)
  weighted_exposure <- pmax(as.numeric(weighted_exposure), 0)
  alpha <- 1 / phi + weighted_count
  beta <- 1 / (phi * prior_mean) + weighted_exposure
  list(mean = alpha / beta, variance = alpha / beta^2,
       alpha = alpha, beta = beta)
}

.dk_p1_stabilized_predict_events <- function(
    x, events, membership, embedding, membership_fit = NULL,
    hard_stratum = NULL, factor_rank = 5L, factor_features = 2000L,
    factor_ridge = 2, min_feature_observed = 20L,
    min_target_observed = 8L, factor_crossfit_folds = 5L,
    factor_crossfit_seed = 1L, bias_kappa = 10,
    support_adaptive_rank = TRUE, predictor_smoothing = 0.25,
    architecture_geometry = NULL) {
  factor_crossfit_folds <- as.integer(factor_crossfit_folds)[1L]
  factor_crossfit_seed <- as.integer(factor_crossfit_seed)[1L]
  if (!is.finite(factor_crossfit_folds) || factor_crossfit_folds < 1L) {
    stop("factor_crossfit_folds must be >= 1", call. = FALSE)
  }
  if (!is.finite(factor_crossfit_seed)) {
    stop("factor_crossfit_seed must be finite", call. = FALSE)
  }
  if (!is.numeric(factor_ridge) || length(factor_ridge) != 1L ||
      !is.finite(factor_ridge) || factor_ridge < 0) {
    stop("factor_ridge must be >= 0", call. = FALSE)
  }
  if (!is.numeric(bias_kappa) || length(bias_kappa) != 1L ||
      is.na(bias_kappa) || bias_kappa < 0) {
    stop("bias_kappa must be a non-negative number or Inf", call. = FALSE)
  }
  if (!is.numeric(predictor_smoothing) || length(predictor_smoothing) != 1L ||
      !is.finite(predictor_smoothing) || predictor_smoothing < 0 ||
      predictor_smoothing > 1) {
    stop("predictor_smoothing must be in [0,1]", call. = FALSE)
  }
  if (!is.logical(support_adaptive_rank) ||
      length(support_adaptive_rank) != 1L || is.na(support_adaptive_rank)) {
    stop("support_adaptive_rank must be TRUE or FALSE", call. = FALSE)
  }
  nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  z <- .dk_align_embedding(embedding, nm$cells)
  hard_stratum <- as.character(.dk_resolve_local_stratum(
    nm$cells, membership, membership_fit, hard_stratum
  ))
  E <- nrow(events)
  out <- list(
    prediction = rep(NA_real_, E),
    factor_prediction = rep(NA_real_, E),
    prediction_sd = rep(NA_real_, E),
    predictability = rep(NA_real_, E),
    shrinkage = numeric(E),
    n_observed_gene = integer(E),
    factor_rank = integer(E),
    factor_features = integer(E),
    factor_iterations = integer(E),
    factor_converged = logical(E),
    factor_fold = integer(E),
    bias_calibration = numeric(E),
    recovery_method = rep("unavailable", E),
    target_mode = rep("positive", E),
    cell_prediction = rep(NA_real_, E),
    cell_available = rep(FALSE, E),
    n_donors = integer(E),
    bandwidth = rep(NA_real_, E),
    local_positive_mean = rep(NA_real_, E),
    local_positive_variance = rep(NA_real_, E),
    local_positive_prevalence = rep(NA_real_, E),
    effective_donors = rep(NA_real_, E),
    tree_distance_weighted_mean = rep(NA_real_, E),
    embedding_distance_weighted_mean = rep(NA_real_, E),
    geometry = architecture_geometry
  )
  if (!E) return(out)

  if (is.null(architecture_geometry)) {
    architecture_geometry <- .dk_build_architecture_geometry(
      z, membership, membership_fit, hard_stratum
    )
  }
  out$geometry <- architecture_geometry

  event_block <- architecture_geometry$same$block_id[events$j]
  for (bid in unique(event_block)) {
    qev <- which(event_block == bid)
    cells <- which(architecture_geometry$same$block_id == bid)
    W <- architecture_geometry$same$W[cells, cells, drop = FALSE]
    genes <- unique(events$i[qev])
    gene_fold <- .dk_target_fold(
      genes, factor_crossfit_folds, factor_crossfit_seed
    )
    for (fold in unique(gene_fold)) {
      fold_genes <- genes[gene_fold == fold]
      qfold <- qev[events$i[qev] %in% fold_genes]
      ev_local <- data.frame(
        i = events$i[qfold], j = match(events$j[qfold], cells)
      )
      fit <- .dk_membership_factor_scores(
        x, cells, ev_local, rank = factor_rank,
        feature_max = factor_features,
        min_feature_observed = min_feature_observed,
        predictor_smoothing = predictor_smoothing,
        smoothing_W = W
      )
      scores <- if (is.null(fit)) {
        matrix(numeric(), length(cells), 0L)
      } else {
        fit$scores
      }
      for (start in seq.int(1L, length(fold_genes), by = 256L)) {
        gset <- fold_genes[start:min(length(fold_genes), start + 255L)]
        y <- as.matrix(x[gset, cells, drop = FALSE])
        state <- dk_batch_p1_state_cpp(
          y, scores, factor_ridge, as.integer(min_target_observed),
          support_adaptive_rank, bias_kappa
        )
        q <- qfold[events$i[qfold] %in% gset]
        gi <- match(events$i[q], gset)
        cj <- match(events$j[q], cells)
        ij <- cbind(gi, cj)
        method_code <- state$method[gi]
        out$prediction[q] <- state$prediction[ij]
        out$prediction_sd[q] <- state$prediction_sd[ij]
        out$shrinkage[q] <- state$shrinkage[gi]
        out$n_observed_gene[q] <- state$n_observed[gi]
        out$n_donors[q] <- state$n_observed[gi]
        out$factor_rank[q] <- state$factor_rank_used[gi]
        out$factor_fold[q] <- fold
        out$bias_calibration[q] <- state$bias_calibration[gi]
        out$recovery_method[q] <- c(
          "unavailable", "positive_membership_mean", "p1_stabilized_state"
        )[method_code + 1L]
        if (!is.null(fit)) {
          out$factor_features[q] <- fit$n_features
          out$factor_iterations[q] <- fit$iterations
          out$factor_converged[q] <- fit$converged
        }
      }
    }
  }
  out
}

.dk_architecture_p1_family <- function(x, counts, events, membership,
                                        geometry, factor_rank = 5L,
                                        factor_features = 2000L,
                                        factor_ridge = 2,
                                        min_feature_observed = 20L,
                                        min_target_observed = 8L,
                                        folds = 5L, seed = 1L,
                                        bias_kappa = 10,
                                        residual_rho = 0.05,
                                        support_kappa = 5,
                                        predictor_smoothing = 0.25,
                                        phi_kappa = 10,
                                        mu_floor = 1e-8,
                                        phi_floor = 1e-4,
                                        scale_factor = 1e4) {
  E <- nrow(events)
  methods <- c("current_p1_crossfit", "p1_local_residual",
               "p1_tree_loo_blend", "p1_stabilized_state",
               "count_eb_p1_local")
  pred <- sd <- matrix(NA_real_, E, length(methods), dimnames = list(NULL, methods))
  diag <- data.frame(event = seq_len(E), effective_donors = NA_real_,
                     residual_alpha = NA_real_, blend_alpha = NA_real_,
                     phi = NA_real_, factor_fold = NA_integer_)
  lib <- as.numeric(Matrix::colSums(counts)); mean_lib <- mean(lib); sf <- lib / mean_lib
  phi_stratum <- .dk_stratum_phi_fallback(counts, sf, geometry$hard_stratum, phi_floor)
  event_block <- geometry$same$block_id[events$j]
  for (bid in unique(event_block)) {
    qev <- which(event_block == bid); cells <- which(geometry$same$block_id == bid)
    W <- geometry$same$W[cells, cells, drop = FALSE]
    W2 <- geometry$same$W2[cells, cells, drop = FALSE]
    genes <- unique(events$i[qev]); gene_fold <- .dk_target_fold(genes, folds, seed)
    for (fold in unique(gene_fold)) {
      fold_genes <- genes[gene_fold == fold]
      qfold <- qev[events$i[qev] %in% fold_genes]
      ev_local <- data.frame(i = events$i[qfold], j = match(events$j[qfold], cells))
      fit <- .dk_membership_factor_scores(
        x, cells, ev_local, rank = factor_rank, feature_max = factor_features,
        min_feature_observed = min_feature_observed
      )
      fit_s <- .dk_membership_factor_scores(
        x, cells, ev_local, rank = factor_rank, feature_max = factor_features,
        min_feature_observed = min_feature_observed,
        predictor_smoothing = predictor_smoothing, smoothing_W = W
      )
      scores <- if (is.null(fit)) matrix(numeric(), length(cells), 0L) else fit$scores
      scores_s <- if (is.null(fit_s)) matrix(numeric(), length(cells), 0L) else fit_s$scores
      for (start in seq.int(1L, length(fold_genes), by = 256L)) {
        gset <- fold_genes[start:min(length(fold_genes), start + 255L)]
        y <- as.matrix(x[gset, cells, drop = FALSE])
        yc <- as.matrix(counts[gset, cells, drop = FALSE])
        st <- dk_batch_p1_state_cpp(y, scores, factor_ridge,
                                    as.integer(min_target_observed), TRUE, bias_kappa)
        st_s <- dk_batch_p1_state_cpp(y, scores_s, factor_ridge,
                                      as.integer(min_target_observed), TRUE, bias_kappa)
        local <- .dk_local_batch_moments(y, W, W2)
        positive <- y > 0 & is.finite(st$loo)
        error <- y - st$loo; error[!positive] <- 0
        Wt <- Matrix::t(W); W2t <- Matrix::t(W2)
        den <- as.matrix(positive %*% Wt)
        den2 <- as.matrix(positive %*% W2t)
        delta <- as.matrix(error %*% Wt) / pmax(den, .Machine$double.eps)
        delta[den <= 0] <- NA_real_
        neff <- matrix(0, nrow(y), ncol(y)); okn <- den2 > 0
        neff[okn] <- den[okn]^2 / den2[okn]
        alpha_res <- alpha_blend <- numeric(nrow(y))
        sigma_res <- sigma_blend <- rep(NA_real_, nrow(y))
        for (gg in seq_len(nrow(y))) {
          donor <- which(positive[gg, ] & is.finite(delta[gg, ]))
          if (length(donor)) {
            e <- error[gg, donor]; dl <- delta[gg, donor]
            ar <- sum(dl * e) / (sum(dl^2) + residual_rho * sum(e^2) + 1e-12)
            alpha_res[gg] <- max(0, min(1, ar))
            si <- neff[gg, donor] / (neff[gg, donor] + support_kappa)
            final_loo <- st$loo[gg, donor] + alpha_res[gg] * si * dl
            sigma_res[gg] <- mean((y[gg, donor] - final_loo)^2)
          }
          donor_b <- which(positive[gg, ] & is.finite(local$mean[gg, ]))
          if (length(donor_b)) {
            d <- local$mean[gg, donor_b] - st$loo[gg, donor_b]
            e <- y[gg, donor_b] - st$loo[gg, donor_b]
            ab <- sum(d * e) / (sum(d^2) + residual_rho * sum(e^2) + 1e-12)
            alpha_blend[gg] <- max(0, min(1, ab))
            si <- local$effective_n[gg, donor_b] /
              (local$effective_n[gg, donor_b] + support_kappa)
            final_loo <- st$loo[gg, donor_b] + alpha_blend[gg] * si * d
            sigma_blend[gg] <- mean((y[gg, donor_b] - final_loo)^2)
          }
        }
        q <- qfold[events$i[qfold] %in% gset]
        gi <- match(events$i[q], gset); cj <- match(events$j[q], cells)
        ij <- cbind(gi, cj)
        base <- st$prediction[ij]; base_sd <- st$prediction_sd[ij]
        pred[q, "current_p1_crossfit"] <- base
        sd[q, "current_p1_crossfit"] <- base_sd
        sr <- neff[ij] / (neff[ij] + support_kappa)
        pr <- base + alpha_res[gi] * sr * delta[ij]
        pr[!is.finite(pr)] <- base[!is.finite(pr)]
        pred[q, "p1_local_residual"] <- pmax(pr, 0)
        vr <- sigma_res[gi] * (1 + 1 / pmax(neff[ij], 1))
        vr[!is.finite(vr) | vr < 0] <- base_sd[!is.finite(vr) | vr < 0]^2
        sd[q, "p1_local_residual"] <- sqrt(pmax(vr, 0))
        sb <- local$effective_n[ij] / (local$effective_n[ij] + support_kappa)
        pb <- base + alpha_blend[gi] * sb * (local$mean[ij] - base)
        pb[!is.finite(pb)] <- base[!is.finite(pb)]
        pred[q, "p1_tree_loo_blend"] <- pmax(pb, 0)
        vb <- sigma_blend[gi] * (1 + 1 / pmax(local$effective_n[ij], 1))
        vb[!is.finite(vb) | vb < 0] <- base_sd[!is.finite(vb) | vb < 0]^2
        sd[q, "p1_tree_loo_blend"] <- sqrt(pmax(vb, 0))
        pred[q, "p1_stabilized_state"] <- st_s$prediction[ij]
        sd[q, "p1_stabilized_state"] <- st_s$prediction_sd[ij]

        u <- sweep(yc, 2L, sf[cells], "/")
        mu_loo <- pmax(expm1(st$loo) / scale_factor * mean_lib, mu_floor)
        phi_hat <- rep(phi_floor, nrow(y))
        for (gg in seq_len(nrow(y))) {
          donor <- which(positive[gg, ] & is.finite(mu_loo[gg, ]))
          if (length(donor)) {
            numerator <- max(sum((u[gg, donor] - mu_loo[gg, donor])^2 -
                                   mu_loo[gg, donor] / sf[cells[donor]]), 0)
            phi_hat[gg] <- max(phi_floor,
              numerator / (sum(mu_loo[gg, donor]^2) + 1e-12))
          }
        }
        stratum_name <- geometry$hard_stratum[cells[1L]]
        phi_s <- phi_stratum[gset, match(stratum_name, colnames(phi_stratum))]
        npos <- rowSums(positive)
        rr <- npos / (npos + phi_kappa)
        phi <- pmax(rr * phi_hat + (1 - rr) * phi_s, phi_floor)
        mu_prior <- pmax(expm1(base) / scale_factor * mean_lib, mu_floor)
        count_num <- as.matrix(yc %*% Wt)[ij]
        sf_num <- as.matrix((yc > 0) %*% Matrix::Diagonal(x = sf[cells]) %*% Wt)[ij]
        denq <- local$positive_weight[ij]; neq <- local$effective_n[ij]
        count_evidence <- ifelse(denq > 0, neq * count_num / denq, 0)
        sf_evidence <- ifelse(denq > 0, neq * sf_num / denq, 0)
        posterior <- .dk_gamma_power_posterior(
          mu_prior, phi[gi], count_evidence, sf_evidence,
          mu_floor = mu_floor, phi_floor = phi_floor
        )
        lambda_post <- posterior$mean
        var_lambda <- posterior$variance
        pe <- log1p(scale_factor * lambda_post / mean_lib)
        deriv <- (scale_factor / mean_lib) /
          (1 + scale_factor * lambda_post / mean_lib)
        pred[q, "count_eb_p1_local"] <- pmax(pe, 0)
        sd[q, "count_eb_p1_local"] <- sqrt(pmax(deriv^2 * var_lambda, 0))
        diag$effective_donors[q] <- neff[ij]
        diag$residual_alpha[q] <- alpha_res[gi]
        diag$blend_alpha[q] <- alpha_blend[gi]
        diag$phi[q] <- phi[gi]
        diag$factor_fold[q] <- fold
      }
    }
  }
  list(prediction = pred, prediction_sd = sd, diagnostics = diag)
}

.dk_hierarchical_partial_pool <- function(x, events, membership, hard_stratum,
                                           geometry, p1_prediction, p1_sd,
                                           kappa_same = 5,
                                           kappa_sibling = 10) {
  mem <- .dk_positive_group_moments(x, membership)
  str <- .dk_positive_group_moments(x, hard_stratum)
  E <- nrow(events); mu1 <- var1 <- rep(NA_real_, E); n1 <- numeric(E)
  Ws <- geometry$W_sibling; W2s <- geometry$W2_sibling
  for (start in seq.int(1L, length(unique(events$i)), by = 256L)) {
    ug <- unique(events$i); gset <- ug[start:min(length(ug), start + 255L)]
    q <- which(events$i %in% gset)
    y <- x[gset, , drop = FALSE]
    lm <- .dk_local_batch_moments(y, Ws, W2s)
    ij <- cbind(match(events$i[q], gset), events$j[q])
    mu1[q] <- lm$mean[ij]; var1[q] <- lm$variance[ij]; n1[q] <- lm$effective_n[ij]
  }
  mi <- match(membership[events$j], mem$levels)
  si <- match(hard_stratum[events$j], str$levels)
  idx0 <- cbind(events$i, mi); idx2 <- cbind(events$i, si)
  mu0 <- mem$mean[idx0]; var0 <- mem$variance[idx0]; n0 <- mem$n[idx0]
  mu2 <- str$mean[idx2]; var2 <- str$variance[idx2]; n2 <- str$n[idx2]
  mu0[!is.finite(mu0)] <- mu2[!is.finite(mu0)]
  var0[!is.finite(var0)] <- var2[!is.finite(var0)]
  mu1[!is.finite(mu1)] <- mu2[!is.finite(mu1)]
  var1[!is.finite(var1)] <- var2[!is.finite(var1)]
  r0 <- n0 / (n0 + kappa_same); r1 <- n1 / (n1 + kappa_sibling)
  r0[!is.finite(r0)] <- 0; r1[!is.finite(r1)] <- 0
  mu_h <- r0 * mu0 + (1 - r0) * (r1 * mu1 + (1 - r1) * mu2)
  # The document's extra outer (1-r0) would square the low-support shrinkage.
  # Replacing the P1 amplitude anchor once is the coherent partial-pooling form.
  prediction <- pmax(p1_prediction + (mu_h - mu0), 0)
  anchor_var <- r0^2 * var0 / pmax(n0, 1) +
    (1 - r0)^2 * (r1^2 * var1 / pmax(n1, 1) +
                    (1 - r1)^2 * var2 / pmax(n2, 1))
  anchor_var[!is.finite(anchor_var) | anchor_var < 0] <- 0
  prediction_sd <- sqrt(pmax(p1_sd^2 + anchor_var, 0))
  list(prediction = prediction, prediction_sd = prediction_sd,
       r_same = r0, r_sibling = r1, sibling_effective_n = n1)
}

.dk_local_simplex <- function(x, events, membership, geometry, folds = 5L,
                              seed = 1L, reliable_features = 200L,
                              lambda = 1, prune_weight = 0.01,
                              min_feature_observed = 20L) {
  # A single target-safe state per cell is sufficient: every masked target in
  # that cell is zero and is therefore excluded by the query-observed feature
  # gate.  This is stricter than fold-wise exclusion and implements the design's
  # explicit "one solve per cell, shared by all dropout genes" requirement.
  K <- ncol(geometry$candidate_index)
  cell_weights <- matrix(0, ncol(x), K)
  n_features_cell <- integer(ncol(x)); solved_cell <- logical(ncol(x))
  query_cells <- sort(unique(events$j))
  for (m in unique(membership[query_cells])) {
    qcells <- query_cells[membership[query_cells] == m]
    cand_global <- geometry$candidate_index[qcells, , drop = FALSE]
    prior <- geometry$candidate_prior[qcells, , drop = FALSE]
    fit_cells <- sort(unique(c(qcells, cand_global[cand_global > 0])))
    cand_local <- matrix(match(cand_global, fit_cells), nrow(cand_global),
                         ncol(cand_global))
    cand_local[is.na(cand_local)] <- 0L
    query_local <- match(qcells, fit_cells)
    y <- x[, fit_cells, drop = FALSE]
    nobs <- as.numeric(.dk_row_sums(y > 0)); s1 <- as.numeric(.dk_row_sums(y))
    s2 <- as.numeric(.dk_row_sums(y * y)); vv <- pmax(
      s2 / pmax(nobs, 1) - (s1 / pmax(nobs, 1))^2, 0
    )
    good <- which(nobs >= min_feature_observed & is.finite(vv) & vv > 1e-10)
    if (length(good) < 2L) next
    # Five times the requested maximum gives sparse query cells enough observed
    # candidates while keeping each local dense working block bounded.
    pool_n <- min(length(good), max(as.integer(reliable_features) * 5L,
                                    as.integer(reliable_features) + 100L))
    good <- good[head(order(vv[good], decreasing = TRUE), pool_n)]
    raw_features <- as.matrix(y[good, , drop = FALSE])
    query_observed <- raw_features[, query_local, drop = FALSE] > 0
    features <- t(scale(t(raw_features)))
    features[!is.finite(features)] <- 0
    fit <- dk_simplex_cell_weights_cpp(
      features, query_observed, as.integer(query_local), cand_local, prior,
      as.integer(reliable_features), lambda, prune_weight, 1e-10
    )
    cell_weights[qcells, ] <- fit$weights
    n_features_cell[qcells] <- fit$n_features
    solved_cell[qcells] <- fit$solved
  }
  sparse_x <- if (inherits(x, "dgCMatrix")) x else methods::as(x, "dgCMatrix")
  ans <- dk_simplex_predict_events_cpp(
    sparse_x, as.integer(events$i), as.integer(events$j),
    geometry$candidate_index, cell_weights, 1e-8
  )
  c(ans, list(cell_weights = cell_weights, n_features_cell = n_features_cell,
              solved_cell = solved_cell))
}

#' Evaluate target-safe recovery architectures on a supplied dropout mask
#'
#' This expert API predicts only coordinates in `mask`. It jointly builds the
#' cross-fitted P1 state and reuses it across local-residual, LOO blend,
#' stabilized-state, hierarchy, simplex, fast, and count empirical-Bayes routes.
#' Raw UMI `counts` are required for the count route. Observed non-target entries
#' are never modified.
#'
#' @export
recovery_architecture_prediction <- function(
    x, counts, membership, embedding, mask, membership_fit = NULL,
    hard_stratum = NULL, methods = .dk_architecture_methods(),
    factor_rank = 5L, factor_features = 2000L, factor_ridge = 2,
    min_feature_observed = 20L, min_target_observed = 8L,
    factor_crossfit_folds = 5L, factor_crossfit_seed = 1L,
    bias_kappa = 10, architecture_geometry = NULL,
    return_details = TRUE) {
  x <- .dk_validate_expression(x); counts <- .dk_validate_expression(counts)
  if (!identical(dim(x), dim(counts))) stop("x and counts dimensions differ", call. = FALSE)
  nm <- .dk_names(x); membership <- .dk_align_membership(membership, nm$cells)
  z <- .dk_align_embedding(embedding, nm$cells)
  if (inherits(membership_fit, "DropoutKillerMembership")) {
    membership_fit$membership <- membership
  }
  hard_stratum <- as.character(.dk_resolve_local_stratum(
    nm$cells, membership, membership_fit, hard_stratum
  ))
  methods <- match.arg(methods, .dk_architecture_methods(), several.ok = TRUE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) {
    if (!all(as.vector(x[cbind(events$i, events$j)] == 0)))
      stop("mask contains observed non-zero entries", call. = FALSE)
    events$membership <- membership[events$j]
  } else events$membership <- integer()
  if (is.null(architecture_geometry)) architecture_geometry <-
    .dk_build_architecture_geometry(z, membership, membership_fit, hard_stratum)
  fam <- .dk_architecture_p1_family(
    x, counts, events, membership, architecture_geometry,
    factor_rank, factor_features, factor_ridge,
    min_feature_observed, min_target_observed,
    factor_crossfit_folds, factor_crossfit_seed, bias_kappa
  )
  pred <- fam$prediction; sd <- fam$prediction_sd
  if ("fast_p1_lite" %in% methods) {
    fast <- .dk_masked_factor_predict_events(
      x, membership, events, factor_rank, factor_features, factor_ridge,
      min_feature_observed, min_target_observed, NULL, "positive",
      1L, factor_crossfit_seed, TRUE, bias_kappa
    )
    pred <- cbind(pred, fast_p1_lite = fast$prediction)
    sd <- cbind(sd, fast_p1_lite = fast$prediction_sd)
  }
  if ("hierarchical_partial_pool" %in% methods) {
    h <- .dk_hierarchical_partial_pool(
      x, events, membership, hard_stratum, architecture_geometry,
      pred[, "current_p1_crossfit"], sd[, "current_p1_crossfit"]
    )
    pred <- cbind(pred, hierarchical_partial_pool = h$prediction)
    sd <- cbind(sd, hierarchical_partial_pool = h$prediction_sd)
  } else h <- NULL
  if ("local_simplex" %in% methods) {
    sx <- .dk_local_simplex(
      x, events, membership, architecture_geometry,
      factor_crossfit_folds, factor_crossfit_seed
    )
    pred <- cbind(pred, local_simplex = sx$prediction)
    sd <- cbind(sd, local_simplex = sx$prediction_sd)
  } else sx <- NULL
  pred <- pred[, intersect(methods, colnames(pred)), drop = FALSE]
  sd <- sd[, colnames(pred), drop = FALSE]
  result <- list(events = events, prediction = pred, prediction_sd = sd,
                 diagnostics = fam$diagnostics, hierarchy = h,
                 simplex = sx, architecture_geometry = architecture_geometry)
  if (return_details) return(result)
  result$prediction
}
