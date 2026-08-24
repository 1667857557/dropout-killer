test_that("masked factor uses unmasked zeros and falls back to the membership mean", {
  x <- matrix(c(0, 0, 0, 4, 4,
                1, 2, 3, 4, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- paste0("c", 1:5)
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  ev <- masked_factor_prediction(x, rep(1, 5), mask, factor_rank = 1,
                                 factor_features = 2, min_feature_observed = 3,
                                 min_target_observed = 10, return_events = TRUE)
  expect_equal(ev$prediction[1], mean(c(0, 0, 4, 4)))
  expect_equal(ev$predictability[1], 0)
  expect_equal(ev$recovery_method[1], "membership_mean")
  expect_true(is.finite(ev$prediction_sd[1]))
})

test_that("masked factor learns held-out coexpression without using the target value", {
  set.seed(41)
  n <- 80
  state <- seq(-2, 2, length.out = n)
  x <- rbind(
    g1 = pmax(0, 5 + 2.5 * state + rnorm(n, sd = 0.15)),
    g2 = pmax(0, 4 + 2.0 * state + rnorm(n, sd = 0.15)),
    g3 = pmax(0, 3 - 1.5 * state + rnorm(n, sd = 0.15)),
    g4 = pmax(0, 2 + 0.8 * state + rnorm(n, sd = 0.20)),
    g5 = pmax(0, 6 - 0.5 * state + rnorm(n, sd = 0.20))
  )
  colnames(x) <- paste0("c", seq_len(n))
  hold <- seq(5, 75, by = 10)
  truth <- x[1, hold]
  xm <- x; xm[1, hold] <- 0
  mask <- Matrix::sparseMatrix(i = rep(1, length(hold)), j = hold, x = TRUE, dims = dim(xm))
  ev <- masked_factor_prediction(xm, rep(1, n), mask, factor_rank = 2,
                                 factor_features = 5, min_feature_observed = 20,
                                 min_target_observed = 20, return_events = TRUE)
  baseline <- mean(xm[1, -hold])
  expect_lt(sqrt(mean((ev$prediction - truth)^2)), sqrt(mean((baseline - truth)^2)))
  expect_gt(mean(ev$predictability), 0)
  expect_true(all(is.finite(ev$prediction_sd)))
  expect_true(all(ev$factor_iterations[ev$factor_rank > 0] == 1L))
  expect_true(all(ev$factor_converged[ev$factor_rank > 0]))
})

test_that("uncertainty-aware draws preserve all observed coordinates", {
  x <- matrix(c(0, 0, 4, 4,
                1, 2, 3, 4), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- paste0("c", 1:4)
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  d <- recover_dropout_expression(x, mask, rep(1, 4), factor_rank = 1,
                                  factor_features = 2, min_feature_observed = 3,
                                  min_target_observed = 10, return_details = TRUE)
  res <- list(expression = d$expression, events = d$events,
              uncertainty_available = TRUE,
              settings = list(recovery_method = "masked_factor"))
  class(res) <- "DropoutKillerResult"
  draws <- sample_dropout_expression(res, n = 3, seed = 9)
  observed <- x != 0
  expect_equal(length(draws), 3)
  for (z in draws) expect_equal(z[observed], x[observed])
})

test_that("uncertainty sampling rejects engines without an uncertainty model", {
  x <- matrix(c(0, 8, 2), nrow = 1, dimnames = list("g1", c("c1", "c2", "c3")))
  z <- matrix(c(0, 0.2, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  d <- recover_dropout_expression(x, mask, rep(1, 3), z, neighbor_k = 2,
                                  return_details = TRUE, recovery_method = "neighbor")
  res <- list(expression = d$expression, events = d$events,
              uncertainty_available = FALSE,
              predictive_variance = NULL,
              settings = list(recovery_method = "neighbor"))
  class(res) <- "DropoutKillerResult"
  expect_error(sample_dropout_expression(res), "unavailable")
})
