test_that("selection can never overwrite observed non-zero values", {
  x <- matrix(c(5, 0, 4, 3, 0, 2), nrow = 2)
  score <- matrix(0.99, nrow(x), ncol(x))
  mask <- select_dropout_mask(score, threshold = 0.95, x = x)
  idx <- which(as.matrix(mask), arr.ind = TRUE)
  expect_true(all(x[idx] == 0))
})

test_that("default finite-sample detector returns sparse zero-only scores", {
  set.seed(3)
  x <- matrix(rexp(240), 12, 20); x[x < 0.5] <- 0
  colnames(x) <- paste0("c", 1:20); rownames(x) <- paste0("g", 1:12)
  s <- local_alra_score(x, rep(1, 20), rank = 2, min_cells = 4)
  expect_true(inherits(s, "sparseMatrix"))
  sm <- Matrix::summary(s)
  if (nrow(sm)) expect_true(all(x[cbind(sm$i, sm$j)] == 0))
  det <- attr(s, "detection")
  expect_equal(det$settings$detection_method, "eb_zero_null")
})

test_that("negative-tail variance shrinkage stays finite with few negative values", {
  yh <- rbind(
    g1 = c(-2, -1, 0.2, 0.3, 0.4, 0.5),
    g2 = c(-0.2, 0.1, 0.2, 0.3, 0.4, 0.5),
    g3 = c(-1.0, -0.8, -0.6, 0.2, 0.3, 0.4)
  )
  z <- DropoutKiller:::.dk_shrunk_negative_null(yh, min_negative = 3,
                                                 variance_prior_df = 10)
  expect_true(all(is.finite(z$sigma)))
  expect_gt(z$sigma[2], 0)
  expect_lt(z$shrinkage_weight[2], z$shrinkage_weight[1])
  expect_true(is.finite(z$prior_sigma))
})

test_that("EB zero-null confidence is one minus gene-wise BH q value", {
  y <- matrix(0, nrow = 2, ncol = 12)
  yh <- rbind(
    c(-1.2, -0.8, -0.5, -0.3, 4.5, 3.8, 0.2, 0.1, 0.3, 0.2, 0.1, 0.2),
    c(-1.0, -0.7, -0.4, -0.2, 2.8, 2.4, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2)
  )
  d <- DropoutKiller:::.dk_detect_eb_zero_null(
    y, yh, seq_len(ncol(y)), 1L, list(rank = 2L), 3L, 10
  )
  expect_true(nrow(d$events) > 0)
  expect_equal(d$events$confidence, 1 - d$events$q_value)
  expect_true(all(d$events$q_value >= d$events$p_value - 1e-12))
  expect_true(all(d$events$null_sigma > 0))
})

test_that("historical ALRA quantile detector remains available", {
  set.seed(4)
  x <- matrix(rexp(240), 12, 20); x[x < 0.5] <- 0
  colnames(x) <- paste0("c", 1:20); rownames(x) <- paste0("g", 1:12)
  d <- local_alra_detect(x, rep(1, 20), rank = 2, quantile_prob = 0.1,
                         min_cells = 4, detection_method = "alra_quantile")
  expect_equal(d$settings$detection_method, "alra_quantile")
  expect_true(all(c("threshold", "confidence") %in% names(d$events)))
})
