test_that("scGACL Gamma-Normal EM matches the released Python implementation", {
  raw <- c(0, 0, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34)
  xlog <- log(1.01 + raw)
  fit <- DropoutKiller:::.dk_scgacl_fit_gene(
    xlog, point = log(1.01)
  )

  # Frozen from Hhjyl/scGACL evaluation/dropout_identify.py using SciPy.
  ref <- c(
    rate = 0.36135584,
    alpha = 0.30264197,
    beta = 0.54655100,
    mu = 1.93044231,
    sigma = 1.00254536
  )
  expect_equal(fit$status, "ok")
  expect_equal(fit$iterations, 2L)
  expect_equal(as.numeric(fit$params), unname(ref), tolerance = 1e-5)

  d0 <- DropoutKiller:::.dk_scgacl_dropout_probability(
    log(1.01), fit$params
  )
  expect_equal(d0, 0.9841226020917497, tolerance = 1e-6)
})

test_that("scGACL invalid-gene rules match released source", {
  # >95% transformed zeros is rejected in get_mix_internal().
  x <- c(rep(log(1.01), 20), log(2.01))
  fit <- DropoutKiller:::.dk_scgacl_fit_gene(x, log(1.01))
  expect_equal(fit$status, "dropout_rate_gt_0.95")
  expect_true(all(is.na(fit$params)))

  counts <- matrix(0, 2, 20)
  counts[1, 1:4] <- c(1, 2, 4, 8)
  counts[2, ] <- 0
  rownames(counts) <- c("expressed", "null")
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  det <- DropoutKiller:::.dk_scgacl_detect(
    counts, group = rep("A", ncol(counts)), dropout_threshold = 0.5
  )
  st <- det$membership_stats
  expect_equal(st$invalid_genes, 1L)
})

test_that("scGACL detector selects only raw zero coordinates within supplied groups", {
  base <- rbind(
    g1 = c(0, 0, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34),
    g2 = c(0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8, 13)
  )
  x <- cbind(base, base)
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  group <- rep(c("A", "B"), each = ncol(base))

  det <- DropoutKiller:::.dk_scgacl_detect(
    x, group = group, dropout_threshold = 0.5
  )
  expect_s3_class(det, "DropoutKillerDetection")
  expect_equal(det$settings$detection_method, "scgacl_gamma_normal")
  expect_equal(det$settings$detection_scope, "group")
  expect_true(nrow(det$events) > 0)
  expect_true(all(x[cbind(det$events$i, det$events$j)] == 0))
  expect_true(all(det$events$confidence >= 0.5))
  expect_setequal(det$events$detection_block, c("A", "B"))

  # Every raw zero for one gene/subpopulation has x=log(1.01), therefore the
  # released mixture assigns an identical posterior to those zero coordinates.
  by_key <- split(det$events$confidence,
                  interaction(det$events$detection_block, det$events$gene, drop = TRUE))
  expect_true(all(vapply(by_key, function(v) length(unique(round(v, 12))) == 1L,
                         logical(1))))
})

test_that("scGACL gene batching and sparse storage do not change detector results", {
  base <- rbind(
    g1 = c(0, 0, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34),
    g2 = c(0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8, 13),
    g3 = c(0, 1, 0, 1, 2, 3, 5, 8, 8, 13, 13, 21),
    g4 = c(0, 0, 2, 2, 3, 3, 4, 5, 7, 9, 12, 18),
    g5 = c(1, 0, 1, 0, 2, 2, 3, 5, 6, 8, 10, 15)
  )
  x <- cbind(base, base)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  group <- rep(c("A", "B"), each = ncol(base))

  dense_one <- DropoutKiller:::.dk_scgacl_detect(
    x, group = group, dropout_threshold = 0.5, gene_batch_size = 1L
  )
  dense_many <- DropoutKiller:::.dk_scgacl_detect(
    x, group = group, dropout_threshold = 0.5, gene_batch_size = 256L
  )
  sparse_many <- DropoutKiller:::.dk_scgacl_detect(
    Matrix::Matrix(x, sparse = TRUE), group = group,
    dropout_threshold = 0.5, gene_batch_size = 256L
  )

  key <- function(d) paste(d$events$i, d$events$j, sep = ":")
  expect_equal(key(dense_one), key(dense_many))
  expect_equal(key(dense_many), key(sparse_many))
  expect_equal(dense_one$events$confidence, dense_many$events$confidence,
               tolerance = 1e-12)
  expect_equal(dense_many$events$confidence, sparse_many$events$confidence,
               tolerance = 1e-12)
  expect_equal(dense_one$events$mixture_rate, dense_many$events$mixture_rate,
               tolerance = 1e-12)
  expect_equal(dense_one$events$gamma_shape, dense_many$events$gamma_shape,
               tolerance = 1e-12)
  expect_equal(dense_one$membership_stats$zero_tests,
               dense_many$membership_stats$zero_tests)
  expect_equal(dense_many$membership_stats$zero_tests,
               sparse_many$membership_stats$zero_tests)
})

test_that("high-level default is scGACL and detector is independent of recovery normalization", {
  base <- rbind(
    g1 = c(0, 0, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34),
    g2 = c(0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8, 13),
    g3 = c(1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 10)
  )
  x <- cbind(base, base)
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  z <- matrix(seq_len(ncol(x) * 2), ncol(x), 2,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  group <- rep(c("A", "B"), each = ncol(base))
  membership <- rep(1:4, each = 6)

  fit1 <- dropout_killer(
    x, z, membership = membership, group = group,
    recovery_method = "neighbor", normalization_scale_factor = 1e4
  )
  fit2 <- dropout_killer(
    x, z, membership = membership, group = group,
    recovery_method = "neighbor", normalization_scale_factor = 5e3
  )

  expect_equal(fit1$settings$detection_method, "scgacl_gamma_normal")
  expect_equal(fit1$settings$detection_scope, "group")
  expect_equal(fit1$settings$detection_input, "raw_counts_log1.01")
  expect_equal(as.matrix(fit1$mask), as.matrix(fit2$mask))
  expect_equal(fit1$detection$events$confidence,
               fit2$detection$events$confidence, tolerance = 0)
})

test_that("default scGACL requires source-faithful subpopulation labels", {
  x <- matrix(c(0, 1, 2, 3, 0, 1, 2, 4), 2, 4)
  rownames(x) <- c("g1", "g2")
  colnames(x) <- paste0("c", 1:4)
  z <- matrix(seq_len(8), 4, 2,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  expect_error(
    dropout_killer(x, z, membership = rep(1, 4), recovery_method = "neighbor"),
    "requires group labels"
  )
})

test_that("scGACL high-level path rejects pre-normalized non-count input", {
  x <- matrix(c(0, 0.2, 1.4, 2.7, 0, 0.4, 1.2, 3.1), 2, 4)
  rownames(x) <- c("g1", "g2")
  colnames(x) <- paste0("c", 1:4)
  z <- matrix(seq_len(8), 4, 2,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  expect_error(
    dropout_killer(
      x, z, membership = rep(1, 4), group = rep("A", 4),
      recovery_method = "neighbor"
    ),
    "requires raw count values"
  )
})
