.dk_scgacl_point <- function() log(1.01)

.dk_scgacl_weight <- function(xdata, params) {
  rate <- params[1L]
  alpha <- params[2L]
  beta <- params[3L]
  mu <- params[4L]
  sigma <- params[5L]
  pz1 <- rate * stats::dgamma(xdata, shape = alpha, rate = beta)
  pz2 <- (1 - rate) * stats::dnorm(xdata, mean = mu, sd = sigma)
  pz <- pz1 / (pz1 + pz2)
  pz[pz1 == 0] <- 0
  list(pz, 1 - pz)
}

.dk_scgacl_update_gamma <- function(xdata, wt) {
  tp_s <- sum(wt)
  tp_t <- sum(wt * xdata)
  tp_u <- sum(wt * log(xdata))
  tp_v <- -tp_u / tp_s - log(tp_s / tp_t)

  if (!is.finite(tp_s) || !is.finite(tp_t) || tp_s <= 0 || tp_t <= 0 ||
      !is.finite(tp_v)) {
    return(c(alpha = NA_real_, beta = NA_real_))
  }

  if (tp_v <= 0) {
    alpha <- 20
  } else {
    alpha0 <- (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / (12 * tp_v)
    if (!is.finite(alpha0)) {
      return(c(alpha = NA_real_, beta = NA_real_))
    }
    if (alpha0 >= 20) {
      alpha <- 20
    } else {
      # scGACL's SciPy implementation solves
      # log(alpha) - digamma(alpha) = tp_v, initialized around alpha0.
      # The equation has one positive root, so uniroot is mathematically
      # equivalent to the released scipy.optimize.root call while avoiding a
      # dependency on Python/SciPy inside the R package.
      f <- function(a) log(a) - digamma(a) - tp_v
      alpha <- tryCatch(
        stats::uniroot(f, interval = c(1e-10, 20), tol = 1e-12)$root,
        error = function(e) NA_real_
      )
    }
  }
  beta <- tp_s / tp_t * alpha
  c(alpha = alpha, beta = beta)
}

.dk_scgacl_density <- function(xdata, params) {
  rate <- params[1L]
  alpha <- params[2L]
  beta <- params[3L]
  mu <- params[4L]
  sigma <- params[5L]
  rate * stats::dgamma(xdata, shape = alpha, rate = beta) +
    (1 - rate) * stats::dnorm(xdata, mean = mu, sd = sigma)
}

.dk_scgacl_fit_gene <- function(xdata, point = .dk_scgacl_point()) {
  # Direct mathematical port of
  # Hhjyl/scGACL/evaluation/dropout_identify.py::get_mix_internal().
  # xdata is already log(1.01 + raw_count).
  rate <- sum(xdata == point) / length(xdata)
  if (rate > 0.95) {
    return(list(params = rep(NA_real_, 5L), iterations = 0L,
                status = "dropout_rate_gt_0.95"))
  }
  if (rate == 0) rate <- 0.01

  alpha <- 1.5
  beta <- 1
  xdata_rm <- xdata[xdata > point]
  if (!length(xdata_rm)) {
    return(list(params = rep(NA_real_, 5L), iterations = 0L,
                status = "no_positive_expression"))
  }
  mu <- mean(xdata_rm)
  # numpy.std() uses ddof=0; do not use stats::sd() here.
  sigma <- sqrt(mean((xdata_rm - mu)^2))
  if (sigma == 0) sigma <- 0.01

  params <- c(rate, alpha, beta, mu, sigma)
  eps <- 10
  iter <- 0L
  loglik_old <- 0

  # Official code stops when squared log10-likelihood change <= 0.5 or after
  # iter > 100. isTRUE() also stops when the comparison becomes NA/NaN, matching
  # Python's false while-condition for a NaN scalar.
  while (isTRUE(eps > 0.5)) {
    wt <- .dk_scgacl_weight(xdata, params)
    tp_sum <- c(sum(wt[[1L]]), sum(wt[[2L]]))
    if (any(!is.finite(tp_sum)) || any(tp_sum <= 0)) {
      params[] <- NA_real_
      break
    }

    rate <- tp_sum[1L] / length(wt[[1L]])
    mu <- sum(wt[[2L]] * xdata) / sum(wt[[2L]])
    sigma <- sqrt(sum(wt[[2L]] * ((xdata - mu)^2) / sum(wt[[2L]])))
    gb <- .dk_scgacl_update_gamma(xdata, wt[[1L]])
    alpha <- gb[["alpha"]]
    beta <- gb[["beta"]]
    params <- c(rate, alpha, beta, mu, sigma)

    if (any(!is.finite(params)) || sigma < 0) {
      params[] <- NA_real_
      break
    }
    new_density <- .dk_scgacl_density(xdata, params)
    loglik <- sum(log10(new_density))
    eps <- (loglik - loglik_old)^2
    loglik_old <- loglik
    iter <- iter + 1L
    if (iter > 100L) break
  }

  status <- if (all(is.finite(params))) "ok" else "numerical_failure"
  list(params = params, iterations = iter, status = status)
}

.dk_scgacl_dropout_probability <- function(xdata, params) {
  rate <- params[1L]
  alpha <- params[2L]
  beta <- params[3L]
  mu <- params[4L]
  sigma <- params[5L]
  gam <- rate * stats::dgamma(xdata, shape = alpha, rate = beta)
  nor <- (1 - rate) * stats::dnorm(xdata, mean = mu, sd = sigma)
  out <- gam / (gam + nor)
  # Released cluster_get_dropout_rate() calls np.nan_to_num(..., nan=0)
  # before returning the posterior matrix.
  out[!is.finite(out)] <- 0
  out
}

.dk_empty_scgacl_events <- function() {
  data.frame(
    i = integer(), j = integer(), gene = character(), cell = character(),
    membership = integer(), detection_block = character(), lowrank = numeric(),
    threshold = numeric(), null_sigma = numeric(), z_score = numeric(),
    p_value = numeric(), q_value = numeric(), confidence = numeric(),
    confidence_fallback = logical(), variance_weight = numeric(),
    alra_margin = numeric(), mixture_rate = numeric(), gamma_shape = numeric(),
    gamma_rate = numeric(), normal_mean = numeric(), normal_sd = numeric(),
    stringsAsFactors = FALSE
  )
}

.dk_scgacl_detect <- function(x, group, dropout_threshold = 0.5,
                              point = .dk_scgacl_point(),
                              gene_batch_size = 256L) {
  x <- .dk_validate_expression(x)
  nm <- .dk_names(x)
  if (is.null(group)) {
    stop(
      "detection_method='scgacl_gamma_normal' requires group labels. ",
      "The official scGACL configuration uses cell-type subpopulations ",
      "(IDENTIFY_USE_CELLTYPE=TRUE); supply those labels through group, or ",
      "select detection_method='alra_global'.",
      call. = FALSE
    )
  }
  grp <- as.character(.dk_align_vector(group, nm$cells, "group", allow_null = FALSE))
  if (anyNA(grp) || any(!nzchar(grp))) {
    stop("group contains missing or empty scGACL subpopulation labels", call. = FALSE)
  }
  if (!is.numeric(dropout_threshold) || length(dropout_threshold) != 1L ||
      !is.finite(dropout_threshold) || dropout_threshold < 0 || dropout_threshold > 1) {
    stop("scgacl_dropout_threshold must be a finite value in [0,1]", call. = FALSE)
  }
  if (!is.numeric(point) || length(point) != 1L || !is.finite(point) || point <= 0) {
    stop("scGACL point must be a positive finite scalar", call. = FALSE)
  }
  gene_batch_size <- as.integer(gene_batch_size)
  if (length(gene_batch_size) != 1L || is.na(gene_batch_size) || gene_batch_size < 1L) {
    stop("gene_batch_size must be a positive integer", call. = FALSE)
  }

  lev <- unique(grp)
  events <- list()
  stats_out <- vector("list", length(lev))
  e <- 0L

  for (ii in seq_along(lev)) {
    label <- lev[ii]
    ids <- which(grp == label)
    n_genes <- nrow(x)
    block_events <- 0L
    invalid_genes <- 0L
    fitted_genes <- 0L
    zero_tests <- 0L

    # The released Python implementation materializes one cell-subpopulation
    # matrix and applies log(1.01 + raw_count). Here genes are processed in
    # bounded batches. Every gene vector passed to EM is byte-for-byte the same
    # mathematical vector; batching changes only memory use, not the model.
    starts <- seq.int(1L, n_genes, by = gene_batch_size)
    for (b0 in starts) {
      gb <- b0:min(n_genes, b0 + gene_batch_size - 1L)
      raw_block <- as.matrix(x[gb, ids, drop = FALSE])
      block <- log(1.01 + raw_block)
      zero_tests <- zero_tests + sum(block == point)

      # Official get_mix_parameters() marks genes invalid when the transformed
      # mean lies within 1e-2 of point before attempting EM.
      genes_expr <- abs(rowMeans(block) - point)
      null_gene <- genes_expr < 1e-2

      for (gg in seq_along(gb)) {
        g <- gb[gg]
        if (null_gene[gg]) {
          invalid_genes <- invalid_genes + 1L
          next
        }
        xdata <- block[gg, ]
        fit <- .dk_scgacl_fit_gene(xdata, point = point)
        if (!identical(fit$status, "ok")) {
          invalid_genes <- invalid_genes + 1L
          next
        }
        fitted_genes <- fitted_genes + 1L
        pars <- fit$params

        # Only observed zeros can become DropoutKiller events. Under the scGACL
        # transform every raw zero has exactly x=point. Thus one posterior
        # evaluation per gene/subpopulation is algebraically identical to
        # constructing the full dense posterior matrix and then selecting zeros.
        zero_local <- which(xdata == point)
        if (!length(zero_local)) next
        d0 <- .dk_scgacl_dropout_probability(point, pars)

        # The paper states d_ij > rho. The released main.py preserves generated
        # values unless predict_drop < threshold, so equality is selected in the
        # actual code. >= therefore reproduces the released implementation.
        if (!is.finite(d0) || d0 < dropout_threshold) next

        jj <- ids[zero_local]
        nadd <- length(jj)
        e <- e + 1L
        block_events <- block_events + nadd
        events[[e]] <- data.frame(
          i = rep.int(g, nadd),
          j = jj,
          gene = rep.int(nm$genes[g], nadd),
          cell = nm$cells[jj],
          membership = rep.int(NA_integer_, nadd),
          detection_block = rep.int(label, nadd),
          lowrank = rep.int(NA_real_, nadd),
          threshold = rep.int(dropout_threshold, nadd),
          null_sigma = rep.int(NA_real_, nadd),
          z_score = rep.int(NA_real_, nadd),
          p_value = rep.int(NA_real_, nadd),
          q_value = rep.int(NA_real_, nadd),
          confidence = rep.int(d0, nadd),
          confidence_fallback = rep.int(FALSE, nadd),
          variance_weight = rep.int(NA_real_, nadd),
          alra_margin = rep.int(NA_real_, nadd),
          mixture_rate = rep.int(pars[1L], nadd),
          gamma_shape = rep.int(pars[2L], nadd),
          gamma_rate = rep.int(pars[3L], nadd),
          normal_mean = rep.int(pars[4L], nadd),
          normal_sd = rep.int(pars[5L], nadd),
          stringsAsFactors = FALSE
        )
      }
    }

    stats_out[[ii]] <- data.frame(
      membership = ii,
      detection_block = label,
      n_cells = length(ids),
      rank = NA_integer_,
      status = "ok",
      alra_candidates = block_events,
      scgacl_candidates = block_events,
      zero_tests = zero_tests,
      prior_sigma = NA_real_,
      invalid_genes = invalid_genes,
      fitted_genes = fitted_genes,
      detection_method = "scgacl_gamma_normal",
      stringsAsFactors = FALSE
    )
  }

  ev <- if (length(events)) do.call(rbind, events) else .dk_empty_scgacl_events()
  st <- if (length(stats_out)) do.call(rbind, stats_out) else data.frame()
  out <- list(
    events = ev,
    membership_stats = st,
    dimensions = dim(x),
    dimnames = list(nm$genes, nm$cells),
    settings = list(
      detection_method = "scgacl_gamma_normal",
      detection_scope = "group",
      dropout_threshold = dropout_threshold,
      point = point,
      transform = "log(1.01 + raw_count)",
      subpopulation_source = "group",
      gene_batch_size = gene_batch_size
    )
  )
  class(out) <- "DropoutKillerDetection"
  out
}
