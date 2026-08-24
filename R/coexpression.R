.dk_membership_factor_scores <- function(x, cells, events, rank = 5L, feature_max = 2000L,
                                         min_feature_observed = 20L) {
  n <- length(cells); G <- nrow(x)
  rank <- max(1L, as.integer(rank)); feature_max <- max(2L, as.integer(feature_max))
  min_feature_observed <- max(3L, as.integer(min_feature_observed))
  if (n < 3L || G < 2L) return(NULL)
  y <- x[, cells, drop = FALSE]
  miss_n <- integer(G)
  if (nrow(events)) miss_n <- tabulate(events$i, nbins = G)
  nobs <- n - miss_n
  sum1 <- as.numeric(.dk_row_sums(y)); sum2 <- as.numeric(.dk_row_sums(y * y))
  mu <- sum1 / pmax(nobs, 1L)
  vv <- rep(0, G); okv <- nobs > 1L
  vv[okv] <- (sum2[okv] - nobs[okv] * mu[okv]^2) / (nobs[okv] - 1L)
  vv <- pmax(vv, 0)
  good <- which(nobs >= min_feature_observed & is.finite(vv) & vv > 1e-10)
  if (nrow(events)) good <- setdiff(good, unique(events$i))
  if (length(good) < 2L) return(NULL)
  if (length(good) > feature_max) {
    ord <- order(vv[good], decreasing = TRUE)
    good <- good[ord[seq_len(feature_max)]]
  }
  k <- min(rank, n - 1L, length(good) - 1L)
  if (k < 1L) return(NULL)
  z <- as.matrix(y[good, , drop = FALSE])
  z <- sweep(z, 1L, mu[good], "-")
  z <- sweep(z, 1L, sqrt(vv[good]), "/")
  z[!is.finite(z)] <- 0
  min_dim <- min(dim(z))
  k <- min(k, min_dim - 1L)
  if (k < 1L) return(NULL)
  take_scores <- function(sv) {
    if (is.null(sv) || is.null(sv$v) || !length(sv$d) || !ncol(sv$v)) return(NULL)
    m <- min(k, length(sv$d), ncol(sv$v))
    if (m < 1L) return(NULL)
    d <- sv$d[seq_len(m)]; lead <- sv$d[1L]
    if (!is.finite(lead) || lead <= 0) return(NULL)
    tol_d <- sqrt(max(1, lead * lead) * 1e-10)
    keep <- which(is.finite(d) & d > tol_d)
    if (!length(keep)) return(NULL)
    sv$v[, keep, drop = FALSE]
  }
  decomposition <- "irlba"
  scores <- NULL
  if (min_dim <= 64L || k >= floor(min_dim / 2L)) {
    sv <- tryCatch(base::svd(z, nu = 0L, nv = k), error = function(e) NULL)
    scores <- take_scores(sv)
    if (!is.null(scores)) decomposition <- "svd_exact"
  } else {
    v0 <- sin(seq_len(ncol(z)) * 1.618033988749895)
    v0 <- v0 / sqrt(sum(v0 * v0))
    sv <- tryCatch(irlba::irlba(z, nu = 0L, nv = k, v = v0), error = function(e) NULL)
    scores <- take_scores(sv)
  }
  if (is.null(scores)) {
    if (min_dim > 512L) return(NULL)
    sv <- tryCatch(base::svd(z, nu = 0L, nv = k), error = function(e) NULL)
    scores <- take_scores(sv)
    if (is.null(scores)) return(NULL)
    decomposition <- "svd_exact_fallback"
  }
  scores <- base::scale(scores, center = TRUE, scale = TRUE)
  scores[!is.finite(scores)] <- 0
  keep_score <- which(colSums(scores * scores) > 1e-10)
  if (!length(keep_score)) return(NULL)
  scores <- scores[, keep_score, drop = FALSE]
  rownames(scores) <- colnames(x)[cells]
  list(scores = scores, rank = ncol(scores), n_features = length(good),
       iterations = 1L, converged = TRUE, decomposition = decomposition)
}

.dk_ridge_factor_target <- function(xg, observed, scores, query, ridge = 1,
                                    min_target_observed = 20L,
                                    target_mode = c("positive", "all_observed")) {
  target_mode <- match.arg(target_mode)
  donor <- observed
  if (target_mode == "positive") donor <- donor & is.finite(xg) & xg > 0
  nobs <- sum(donor); y <- xg[donor]
  mu <- if (nobs) mean(y) else 0
  total_sse <- if (nobs) sum((y - mu)^2) else 0
  total_var <- if (nobs > 1L) total_sse / (nobs - 1L) else NA_real_
  fallback <- function() {
    if (nobs < 2L || !is.finite(total_var)) {
      return(list(
        prediction = rep(0, length(query)), factor_prediction = rep(NA_real_, length(query)),
        prediction_sd = rep(NA_real_, length(query)), predictability = 0, shrinkage = 0,
        n_observed = nobs, method = "unavailable", effective_df = 0,
        target_mode = target_mode
      ))
    }
    pv <- total_var * (1 + 1 / nobs)
    method <- if (target_mode == "positive") "positive_membership_mean" else "membership_mean"
    list(
      prediction = rep(max(0, mu), length(query)), factor_prediction = rep(NA_real_, length(query)),
      prediction_sd = rep(sqrt(max(0, pv)), length(query)), predictability = 0, shrinkage = 0,
      n_observed = nobs, method = method, effective_df = 1, target_mode = target_mode
    )
  }
  if (!nobs || is.null(scores) || !ncol(scores) ||
      nobs < max(as.integer(min_target_observed), ncol(scores) + 3L) ||
      !is.finite(total_sse) || total_sse <= 1e-12) return(fallback())
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0)
    stop("factor_ridge must be >= 0", call. = FALSE)

  Xo <- cbind(1, scores[donor, , drop = FALSE])
  P <- diag(c(0, rep(ridge, ncol(scores))), nrow = ncol(Xo))
  XtX <- crossprod(Xo); A <- XtX + P
  inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(inv) || any(!is.finite(inv))) return(fallback())

  beta <- as.vector(inv %*% crossprod(Xo, y))
  fitted <- as.vector(Xo %*% beta)
  Hdiag <- rowSums((Xo %*% inv) * Xo)
  loo_denom <- pmax(1 - Hdiag, 1e-6)
  loo_factor <- y - (y - fitted) / loo_denom
  loo_null <- if (nobs > 1L) (sum(y) - y) / (nobs - 1L) else rep(mu, nobs)

  d <- loo_factor - loo_null
  target <- y - loo_null
  den <- sum(d * d)
  q <- if (is.finite(den) && den > 1e-12) sum(d * target) / den else 0
  q <- max(0, min(1, q))
  loo_shrunk <- loo_null + q * d

  null_loo_sse <- sum(target^2)
  model_loo_sse <- sum((y - loo_shrunk)^2)
  predictability <- if (is.finite(null_loo_sse) && null_loo_sse > 0)
    1 - model_loo_sse / null_loo_sse else 0
  predictability <- max(0, min(1, predictability))

  Xq <- cbind(1, scores[query, , drop = FALSE])
  raw <- as.vector(Xq %*% beta)
  pred <- pmax(mu + q * (raw - mu), 0)

  df <- sum(diag(inv %*% XtX))
  df_shrunk <- 1 + q * max(df - 1, 0)

  sigma2 <- mean((y - loo_shrunk)^2)
  if (!is.finite(sigma2) || sigma2 < 0) sigma2 <- total_var
  h <- rowSums((Xq %*% inv) * Xq)
  mean_leverage <- (1 - q)^2 / nobs + q^2 * pmax(h, 0)
  pv <- sigma2 * (1 + mean_leverage)
  pv[!is.finite(pv) | pv < 0] <- total_var * (1 + 1 / nobs)

  method <- if (q > 0) {
    if (target_mode == "positive") "masked_factor_positive" else "masked_factor"
  } else {
    if (target_mode == "positive") "positive_membership_mean" else "membership_mean"
  }
  list(
    prediction = pred, factor_prediction = raw, prediction_sd = sqrt(pmax(pv, 0)),
    predictability = predictability, shrinkage = q, n_observed = nobs,
    method = method, effective_df = df_shrunk, target_mode = target_mode
  )
}

.dk_masked_factor_predict_events <- function(x, membership, events, factor_rank = 5L,
                                             factor_features = 2000L, factor_ridge = 1,
                                             min_feature_observed = 20L,
                                             min_target_observed = 20L,
                                             cap_quantile = NULL,
                                             target_mode = c("positive", "all_observed")) {
  target_mode <- match.arg(target_mode)
  n_ev <- nrow(events)
  out <- list(
    prediction = rep(NA_real_, n_ev), factor_prediction = rep(NA_real_, n_ev),
    prediction_sd = rep(NA_real_, n_ev), predictability = numeric(n_ev),
    shrinkage = numeric(n_ev), n_observed_gene = integer(n_ev),
    factor_rank = integer(n_ev), factor_features = integer(n_ev),
    factor_iterations = integer(n_ev), factor_converged = logical(n_ev),
    recovery_method = rep("unavailable", n_ev), target_mode = rep(target_mode, n_ev)
  )
  if (!n_ev) return(out)
  if (!is.null(cap_quantile)) {
    if (!is.numeric(cap_quantile) || length(cap_quantile) != 1L ||
        !is.finite(cap_quantile) || cap_quantile <= 0.5 || cap_quantile > 1)
      stop("cap_quantile must be in (0.5,1]", call. = FALSE)
  }
  for (m in unique(events$membership)) {
    cells <- which(membership == m); q <- which(events$membership == m)
    evm <- events[q, , drop = FALSE]
    fit <- .dk_membership_factor_scores(
      x, cells, evm, rank = factor_rank, feature_max = factor_features,
      min_feature_observed = min_feature_observed
    )
    scores <- if (is.null(fit)) NULL else fit$scores
    for (g in unique(evm$i)) {
      qg_local <- which(evm$i == g); qg <- q[qg_local]
      query <- match(evm$j[qg_local], cells)
      observed <- rep(TRUE, length(cells)); observed[query] <- FALSE
      xg <- as.numeric(x[g, cells])
      tg <- .dk_ridge_factor_target(
        xg, observed, scores, query, ridge = factor_ridge,
        min_target_observed = min_target_observed, target_mode = target_mode
      )
      pred <- tg$prediction
      if (!is.null(cap_quantile)) {
        donors <- observed
        if (target_mode == "positive") donors <- donors & xg > 0
        cap_values <- xg[donors & is.finite(xg)]
        if (length(cap_values)) {
          cap <- as.numeric(stats::quantile(cap_values, cap_quantile, names = FALSE, type = 8))
          if (is.finite(cap)) pred <- pmin(pred, cap)
        }
      }
      out$prediction[qg] <- pred
      out$factor_prediction[qg] <- tg$factor_prediction
      out$prediction_sd[qg] <- tg$prediction_sd
      out$predictability[qg] <- tg$predictability
      out$shrinkage[qg] <- tg$shrinkage
      out$n_observed_gene[qg] <- tg$n_observed
      out$recovery_method[qg] <- tg$method
      out$target_mode[qg] <- tg$target_mode
      if (!is.null(fit)) {
        out$factor_rank[qg] <- fit$rank; out$factor_features[qg] <- fit$n_features
        out$factor_iterations[qg] <- fit$iterations; out$factor_converged[qg] <- fit$converged
      }
    }
  }
  out
}

#' Masked membership-local coexpression prediction
#'
#' Learns a low-dimensional coexpression state inside each membership from
#' high-variance non-target genes. Every gene carrying a recovery event is
#' excluded from factor-state learning, preventing target leakage.
#'
#' By default, once a supplied zero has already been classified as technical
#' dropout, target-gene magnitude is learned only from reliable positive donor
#' values in the same membership. This estimates
#' `E[X_g | X_g > 0, cell state]` rather than the zero-inflated unconditional
#' membership mean. Unmasked zeros are never modified. Set
#' `target_mode = "all_observed"` to reproduce the previous unconditional target.
#'
#' Exact analytic leave-one-out predictions choose shrinkage toward the
#' corresponding membership mean. Predictive standard deviations use LOO
#' residual error so residual biological variability is not discarded.
#'
#' @export
masked_factor_prediction <- function(x, membership, mask, factor_rank = 5L,
                                     factor_features = 2000L, factor_ridge = 1,
                                     min_feature_observed = 20L,
                                     min_target_observed = 20L,
                                     cap_quantile = NULL, return_events = FALSE,
                                     target_mode = c("positive", "all_observed")) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  target_mode <- match.arg(target_mode)
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x))))
    stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) {
    is_zero <- as.vector(x[cbind(events$i, events$j)] == 0)
    if (!all(is_zero)) stop("mask contains observed non-zero entries", call. = FALSE)
    events$membership <- membership[events$j]
  } else events$membership <- integer()
  fit <- .dk_masked_factor_predict_events(
    x, membership, events, factor_rank, factor_features, factor_ridge,
    min_feature_observed, min_target_observed, cap_quantile, target_mode
  )
  if (return_events) {
    events$prediction <- fit$prediction; events$factor_prediction <- fit$factor_prediction
    events$prediction_sd <- fit$prediction_sd; events$predictability <- fit$predictability
    events$shrinkage <- fit$shrinkage; events$n_observed_gene <- fit$n_observed_gene
    events$factor_rank <- fit$factor_rank; events$factor_features <- fit$factor_features
    events$factor_iterations <- fit$factor_iterations; events$factor_converged <- fit$factor_converged
    events$recovery_method <- fit$recovery_method; events$target_mode <- fit$target_mode
    return(events)
  }
  ok <- is.finite(fit$prediction) & fit$prediction > 0
  .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok],
                     nrow(x), ncol(x), list(nm$genes, nm$cells))
}

.dk_positive_gamma_draw <- function(mu, sd) {
  out <- mu
  good <- is.finite(mu) & mu > 0 & is.finite(sd) & sd > 0
  if (any(good)) {
    variance <- sd[good]^2
    shape <- mu[good]^2 / variance
    scale <- variance / mu[good]
    stable <- is.finite(shape) & shape > 0 & is.finite(scale) & scale > 0
    draw <- mu[good]
    if (any(stable)) {
      draw[stable] <- stats::rgamma(sum(stable), shape = shape[stable], scale = scale[stable])
    }
    out[good] <- draw
  }
  pmax(out, 0)
}

#' Draw uncertainty-aware completed expression matrices
#'
#' Draws only previously recovered dropout events from their event-level
#' predictive approximation. Observed entries remain exactly fixed.
#' For the default positive-conditional recovery target, positive draws use a
#' Gamma moment match with the stored predictive mean and variance, so repeated
#' draws preserve those first two moments while remaining non-negative.
#' Results from engines without an uncertainty model are rejected.
#'
#' @export
sample_dropout_expression <- function(result, n = 1L, seed = NULL) {
  if (!inherits(result, "DropoutKillerResult")) stop("result must be a DropoutKillerResult", call. = FALSE)
  if (!isTRUE(result$uncertainty_available)) stop("predictive uncertainty is unavailable for this result; use recovery_method='masked_factor'", call. = FALSE)
  n <- as.integer(n)[1L]
  if (!is.finite(n) || n < 1L) stop("n must be >= 1", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  ev <- result$events
  use <- if (nrow(ev) && all(c("prediction_sd", "recovered", "changed") %in% names(ev))) which(ev$changed) else integer()
  if (length(use) && any(!is.finite(ev$recovered[use]) | !is.finite(ev$prediction_sd[use]) | ev$prediction_sd[use] < 0))
    stop("recovered events do not contain complete predictive uncertainty", call. = FALSE)
  positive_target <- !is.null(result$settings$factor_target) &&
    identical(result$settings$factor_target, "positive")
  one <- function() {
    out <- result$expression
    if (!length(use)) return(out)
    if (positive_target) {
      draw <- .dk_positive_gamma_draw(ev$recovered[use], ev$prediction_sd[use])
    } else {
      draw <- pmax(0, stats::rnorm(length(use), mean = ev$recovered[use], sd = ev$prediction_sd[use]))
    }
    if (inherits(out, "Matrix")) {
      delta <- .dk_sparse_numeric(ev$i[use], ev$j[use], draw - ev$recovered[use],
                                  nrow(out), ncol(out), dimnames(out))
      out <- out + delta
    } else out[cbind(ev$i[use], ev$j[use])] <- draw
    out
  }
  ans <- lapply(seq_len(n), function(i) one())
  if (n == 1L) ans[[1L]] else ans
}
