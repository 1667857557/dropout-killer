test_that("ALRA comparator is one original all-cell block, not cell-class blocks", {
  set.seed(101)
  x <- matrix(rexp(20 * 12), 20, 12)
  x[x < 0.8] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  group <- rep(c("B", "T"), each = 6)

  det <- DropoutKiller:::.dk_global_alra_detect(
    x, group = group, rank = 2, quantile_prob = 0.001,
    min_cells = 4, seed = 101, K = 100, rank_z = 6,
    noise_start = 80, choose_q = 2, svd_q = 2
  )

  expect_s3_class(det, "DropoutKillerDetection")
  expect_equal(det$settings$detection_method, "alra_global")
  expect_equal(det$settings$detection_scope, "all_cells")
  expect_equal(nrow(det$membership_stats), 1L)
  expect_equal(det$membership_stats$detection_block, "all_cells")
  expect_equal(det$membership_stats$n_cells, ncol(x))
})

test_that("original ALRA automatic rank retains source K validity constraints", {
  set.seed(102)
  x <- matrix(rexp(20 * 8), 20, 8)
  x[x < 0.7] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))

  # KlugerLab/ALRA::choose_k() stops when K > min(dim(A_norm)); the old
  # DropoutKiller group-wise detector silently adapted K, which is no longer an
  # ALRA comparator behavior.
  expect_error(
    DropoutKiller:::.dk_original_alra_detect(
      x, rank = "auto", seed = 102, K = 100,
      noise_start = 80, choose_q = 2, svd_q = 2
    ),
    "K must not exceed"
  )
})

test_that("original ALRA zero gate numerically matches released source calculation", {
  set.seed(601)
  raw <- matrix(rpois(40 * 60, lambda = 2), 40, 60)
  raw[1, colSums(raw) == 0] <- 1
  raw[sample(length(raw), 500)] <- 0
  rownames(raw) <- paste0("g", seq_len(nrow(raw)))
  colnames(raw) <- paste0("c", seq_len(ncol(raw)))
  x <- DropoutKiller:::.dk_alra_library_log(raw)
  A <- t(as.matrix(x))

  # Directly reproduce KlugerLab/ALRA::alra() up to the adaptive zero gate.
  # A is cells x genes, exactly the orientation expected by the original code.
  seed <- 602L
  set.seed(seed)
  ref <- rsvd::rsvd(A, k = 3L, q = 2L)
  lr <- ref$u[, 1:3, drop = FALSE] %*%
    diag(ref$d[1:3], nrow = 3L) %*%
    t(ref$v[, 1:3, drop = FALSE])
  tau <- abs(apply(
    lr, 2L, stats::quantile,
    probs = 0.001, names = FALSE
  ))
  ref_pass <- (A == 0) & sweep(lr, 2L, tau, FUN = ">")

  det <- DropoutKiller:::.dk_original_alra_detect(
    x, rank = 3L, quantile_prob = 0.001,
    seed = seed, svd_q = 2L
  )
  got <- matrix(FALSE, nrow(A), ncol(A))
  if (nrow(det$events)) {
    got[cbind(det$events$j, det$events$i)] <- TRUE
  }

  expect_equal(got, ref_pass)
  if (nrow(det$events)) {
    expect_equal(
      det$events$threshold,
      tau[det$events$i],
      tolerance = 1e-12
    )
    expect_equal(
      det$events$lowrank,
      lr[cbind(det$events$j, det$events$i)],
      tolerance = 1e-12
    )
  }
})

test_that("explicit high-level ALRA uses all cells and threshold is not a second gate", {
  set.seed(103)
  x <- matrix(rexp(24 * 12), 24, 12)
  x[x < 0.8] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  z <- matrix(rnorm(12 * 3), 12, 3,
              dimnames = list(colnames(x), paste0("PC", 1:3)))
  group <- rep(c("class1", "class2"), each = 6)
  membership <- rep(seq_len(4), each = 3)

  fit1 <- dropout_killer(
    x, z, membership = membership, group = group,
    detection_method = "alra_global", rank = 2,
    min_cells = 4, threshold = 0.1, seed = 103,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 4, min_target_observed = 4
  )
  fit2 <- dropout_killer(
    x, z, membership = membership, group = group,
    detection_method = "alra_global", rank = 2,
    min_cells = 4, threshold = 0.999, seed = 103,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 4, min_target_observed = 4
  )

  expect_equal(fit1$settings$detection_method, "alra_global")
  expect_equal(fit1$settings$detection_scope, "all_cells")
  expect_equal(fit1$detection$membership_stats$detection_block, "all_cells")
  expect_equal(as.matrix(fit1$mask), as.matrix(fit2$mask))
  if (nrow(fit1$events)) {
    expect_true(all(fit1$events$detection_block == "all_cells"))
    expect_equal(fit1$events$membership, membership[fit1$events$j])
  }
})

test_that("deprecated group-wise ALRA name maps to original global ALRA", {
  set.seed(104)
  x <- matrix(rexp(20 * 12), 20, 12)
  x[x < 0.8] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  z <- matrix(rnorm(12 * 2), 12, 2,
              dimnames = list(colnames(x), c("PC1", "PC2")))

  expect_warning(
    fit <- dropout_killer(
      x, z, membership = rep(c(1, 2), each = 6),
      group = rep(c("A", "B"), each = 6),
      detection_method = "alra_global_by_group", rank = 2,
      seed = 104, recovery_method = "neighbor"
    ),
    "deprecated"
  )
  expect_equal(fit$settings$detection_method, "alra_global")
  expect_equal(fit$settings$detection_scope, "all_cells")
})

test_that("historical EB detector remains explicitly available", {
  set.seed(105)
  x <- matrix(rexp(20 * 12), 20, 12)
  x[x < 0.8] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  z <- matrix(rnorm(12 * 2), 12, 2,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  fit <- dropout_killer(
    x, z, membership = rep(c(1, 2), each = 6),
    detection_method = "eb_zero_null", rank = 2,
    min_cells = 4, threshold = 0.95, seed = 105,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 4, min_target_observed = 4
  )
  expect_equal(fit$settings$detection_method, "eb_zero_null")
  expect_equal(fit$settings$detection_scope, "membership")
})
