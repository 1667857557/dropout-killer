test_that("analytic ridge LOO matches brute-force refitting", {
  set.seed(12)
  n <- 14L; k <- 3L; ridge <- 1.7
  X <- cbind(1, matrix(rnorm(n * k), n, k)); y <- rnorm(n)
  P <- diag(c(0, rep(ridge, k)))
  inv <- solve(crossprod(X) + P)
  fitted <- as.vector(X %*% inv %*% crossprod(X, y))
  h <- rowSums((X %*% inv) * X)
  analytic <- y - (y - fitted) / (1 - h)
  brute <- vapply(seq_len(n), function(i) {
    keep <- setdiff(seq_len(n), i)
    beta <- solve(crossprod(X[keep, , drop = FALSE]) + P,
                  crossprod(X[keep, , drop = FALSE], y[keep]))
    sum(X[i, ] * beta)
  }, numeric(1))
  expect_lt(max(abs(analytic - brute)), 1e-8)
})

test_that("target genes do not enter their factor state", {
  set.seed(13)
  x <- matrix(rnorm(12 * 30), 12, 30)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  events <- data.frame(i = c(1L, 2L), j = c(1L, 2L))
  a <- DropoutKiller:::.dk_membership_factor_scores(
    x, seq_len(ncol(x)), events, rank = 3, min_feature_observed = 3
  )
  x[1:2, ] <- x[1:2, ] + matrix(rnorm(60, 100, 20), 2, 30)
  b <- DropoutKiller:::.dk_membership_factor_scores(
    x, seq_len(ncol(x)), events, rank = 3, min_feature_observed = 3
  )
  expect_equal(a$scores, b$scores, tolerance = 1e-10)
})

test_that("simplex ridge weights satisfy the KKT constraints", {
  set.seed(14)
  f <- matrix(rnorm(8 * 12), 8, 12)
  gram <- crossprod(f) / nrow(f)
  candidates <- matrix(rep(c(2L, 3L, 4L, 5L), 3), 3, 4, byrow = TRUE)
  query <- c(1L, 6L, 7L)
  candidates[2, ] <- c(1L, 2L, 3L, 4L)
  candidates[3, ] <- c(8L, 9L, 10L, 11L)
  prior <- matrix(0.25, 3, 4)
  w <- DropoutKiller:::dk_simplex_weights_cpp(
    gram, query, candidates, prior, 1, 0, 1e-10
  )
  expect_true(all(w >= -1e-10))
  expect_lt(max(abs(rowSums(w) - 1)), 1e-8)
})

test_that("cell-wise simplex excludes every query-masked target in one solve", {
  set.seed(140)
  features <- matrix(rnorm(8 * 7), 8, 7)
  observed <- matrix(TRUE, 8, 2)
  observed[1:3, 1] <- FALSE
  observed[6:8, 2] <- FALSE
  candidates <- rbind(c(2L, 3L, 4L), c(1L, 3L, 4L))
  prior <- matrix(1 / 3, 2, 3)
  fit1 <- DropoutKiller:::dk_simplex_cell_weights_cpp(
    features, observed, c(1L, 2L), candidates, prior, 5L, 1, 0, 1e-10
  )
  changed <- features
  changed[1:3, ] <- matrix(rnorm(3 * 7, 100, 30), 3, 7)
  fit2 <- DropoutKiller:::dk_simplex_cell_weights_cpp(
    changed, observed, c(1L, 2L), candidates, prior, 5L, 1, 0, 1e-10
  )
  expect_equal(fit1$weights[1, ], fit2$weights[1, ], tolerance = 1e-10)
  expect_equal(length(fit1$solved), 2L)
  expect_true(all(rowSums(fit1$weights[fit1$solved, , drop = FALSE]) > 0.999999))
})

test_that("count power posterior has required limiting behavior", {
  prior <- 4
  no_evidence <- DropoutKiller:::.dk_gamma_power_posterior(prior, 0.5, 0, 0)
  expect_equal(no_evidence$mean, prior, tolerance = 1e-12)
  more_count <- DropoutKiller:::.dk_gamma_power_posterior(prior, 0.5, 20, 2)
  expect_gt(more_count$mean, no_evidence$mean)
  concentrated <- DropoutKiller:::.dk_gamma_power_posterior(
    prior, 1e-12, 20, 2, phi_floor = 1e-12
  )
  expect_equal(concentrated$mean, prior, tolerance = 1e-3)
  expect_true(all(more_count$alpha > 0 & more_count$beta > 0 &
                    more_count$variance >= 0))
})

test_that("architecture prediction is event-only and numerically valid", {
  skip_if_not_installed("Matrix")
  set.seed(15)
  counts <- matrix(rpois(24 * 48, 2), 24, 48)
  counts[, colSums(counts) == 0] <- 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  embedding <- matrix(rnorm(48 * 3), 48, 3,
                      dimnames = list(colnames(counts), NULL))
  membership <- rep(seq_len(4), each = 12)
  label <- rep(c("A", "B"), each = 24)
  positive <- which(counts > 0, arr.ind = TRUE)
  take <- positive[sample(seq_len(nrow(positive)), 60), , drop = FALSE]
  masked <- counts; masked[take] <- 0
  sparse_counts <- Matrix::Matrix(masked, sparse = TRUE)
  x <- DropoutKiller:::.dk_alra_library_log(sparse_counts)
  mask <- Matrix::sparseMatrix(i = take[, 1], j = take[, 2], x = TRUE,
                               dims = dim(x), dimnames = dimnames(x)) > 0
  ans <- recovery_architecture_prediction(
    x, sparse_counts, membership, embedding, mask,
    hard_stratum = label, factor_crossfit_folds = 3,
    min_feature_observed = 3, min_target_observed = 3
  )
  expect_equal(nrow(ans$prediction), nrow(take))
  expect_true(all(is.finite(ans$prediction[, "current_p1_crossfit"])))
  expect_true(all(is.finite(ans$prediction[, "count_eb_p1_local"])))
  ok <- is.finite(ans$simplex$sum_weight_error)
  expect_true(!any(ok) || max(ans$simplex$sum_weight_error[ok]) < 1e-8)
  expect_true(!any(ok) || min(ans$simplex$min_weight[ok]) >= -1e-10)
})

test_that("production P1 stabilized dispatcher matches the audited architecture", {
  skip_if_not_installed("Matrix")
  set.seed(151)
  counts <- matrix(rpois(30 * 54, 2), 30, 54)
  counts[, colSums(counts) == 0] <- 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  embedding <- matrix(
    rnorm(54 * 4), 54, 4,
    dimnames = list(colnames(counts), paste0("PC", 1:4))
  )
  membership <- rep(seq_len(6), each = 9)
  label <- rep(c("A", "B"), each = 27)
  positive <- which(counts > 0, arr.ind = TRUE)
  take <- positive[sample(seq_len(nrow(positive)), 80), , drop = FALSE]
  masked <- counts
  masked[take] <- 0
  sparse_counts <- Matrix::Matrix(masked, sparse = TRUE)
  x <- DropoutKiller:::.dk_alra_library_log(sparse_counts)
  mask <- Matrix::sparseMatrix(
    i = take[, 1], j = take[, 2], x = TRUE,
    dims = dim(x), dimnames = dimnames(x)
  ) > 0

  expert <- recovery_architecture_prediction(
    x, sparse_counts, membership, embedding, mask,
    hard_stratum = label, methods = "p1_stabilized_state",
    factor_rank = 3, factor_features = 20, factor_ridge = 2,
    min_feature_observed = 3, min_target_observed = 3,
    factor_crossfit_folds = 3, factor_crossfit_seed = 7,
    bias_kappa = 10
  )
  production <- recover_dropout_expression(
    x, mask, membership, embedding,
    recovery_method = "p1_stabilized_state",
    hard_stratum = label, factor_rank = 3, factor_features = 20,
    factor_ridge = 2, min_feature_observed = 3,
    min_target_observed = 3, factor_crossfit_folds = 3,
    factor_crossfit_seed = 7, support_adaptive_rank = TRUE,
    bias_kappa = 10, predictor_smoothing = 0.25,
    return_details = TRUE
  )

  expect_equal(
    production$events$recovered,
    expert$prediction[, "p1_stabilized_state"],
    tolerance = 1e-10
  )
  expect_equal(
    production$events$prediction_sd,
    expert$prediction_sd[, "p1_stabilized_state"],
    tolerance = 1e-10
  )
  expect_true(all(is.finite(production$events$recovered)))
  keep <- !as.matrix(mask)
  expect_equal(
    as.numeric(production$expression[keep]),
    as.numeric(x[keep]),
    tolerance = 0
  )
})
