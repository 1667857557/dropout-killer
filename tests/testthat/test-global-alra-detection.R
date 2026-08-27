test_that("global ALRA detection uses major cell classes rather than memberships", {
  set.seed(101)
  x <- matrix(rexp(12 * 8), 12, 8)
  x[x < 0.8] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  group <- rep(c("B", "T"), each = 4)
  det <- DropoutKiller:::.dk_global_alra_detect(
    x, group = group, rank = 2, quantile_prob = 0.001,
    min_cells = 4, seed = 101, K = 100, rank_z = 6,
    noise_start = 80, choose_q = 2, svd_q = 2
  )
  expect_s3_class(det, "DropoutKillerDetection")
  expect_equal(det$settings$detection_method, "alra_global_by_group")
  expect_equal(det$settings$detection_scope, "group")
  expect_setequal(det$membership_stats$detection_block, c("B", "T"))
  expect_equal(sort(det$membership_stats$n_cells), c(4, 4))
})

test_that("global ALRA auto rank adapts K to small cell classes", {
  set.seed(102)
  x <- matrix(rexp(20 * 8), 20, 8)
  x[x < 0.7] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  det <- DropoutKiller:::.dk_global_alra_detect(
    x, group = rep("small_class", 8), rank = "auto",
    min_cells = 8, seed = 102, K = 100, noise_start = 80,
    choose_q = 2, svd_q = 2
  )
  expect_equal(det$membership_stats$status, "ok")
  expect_equal(det$membership_stats$K, 7L)
  expect_true(det$membership_stats$rank >= 1L)
  expect_true(det$membership_stats$rank <= 7L)
})

test_that("high-level default is group-global ALRA and threshold is not a second gate", {
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
    rank = 2, min_cells = 4, threshold = 0.1, seed = 103,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 4, min_target_observed = 4
  )
  fit2 <- dropout_killer(
    x, z, membership = membership, group = group,
    rank = 2, min_cells = 4, threshold = 0.999, seed = 103,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 4, min_target_observed = 4
  )
  expect_equal(fit1$settings$detection_method, "alra_global_by_group")
  expect_equal(fit1$settings$detection_scope, "group")
  expect_setequal(fit1$detection$membership_stats$detection_block, unique(group))
  expect_equal(as.matrix(fit1$mask), as.matrix(fit2$mask))
  if (nrow(fit1$events)) {
    expect_true(all(fit1$events$detection_block %in% unique(group)))
    expect_equal(fit1$events$membership, membership[fit1$events$j])
  }
})

test_that("historical EB detector remains explicitly available", {
  set.seed(104)
  x <- matrix(rexp(20 * 12), 20, 12)
  x[x < 0.8] <- 0
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  z <- matrix(rnorm(12 * 2), 12, 2,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  fit <- dropout_killer(
    x, z, membership = rep(c(1, 2), each = 6),
    detection_method = "eb_zero_null", rank = 2,
    min_cells = 4, threshold = 0.95, seed = 104,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 4, min_target_observed = 4
  )
  expect_equal(fit$settings$detection_method, "eb_zero_null")
  expect_equal(fit$settings$detection_scope, "membership")
})
