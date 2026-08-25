test_that("walktrap hierarchy is retained by the original membership builder", {
  set.seed(1)
  z <- cbind(x = c(rnorm(8, -2, 0.15), rnorm(8, 2, 0.15)),
             y = rnorm(16, 0, 0.1))
  rownames(z) <- paste0("c", seq_len(nrow(z)))
  fit <- build_supercell_membership(z, gamma = 4, k_knn = 3,
                                    method = "walktrap", approximate = FALSE)
  expect_s3_class(fit, "DropoutKillerMembership")
  expect_true(length(fit$hierarchies) == 1L)
  expect_false(is.null(fit$hierarchies[[1L]]))
  expect_false(is.null(fit$tree_indices[[1L]]))
  ti <- fit$tree_indices[[1L]]
  d12 <- DropoutKiller:::.dk_tree_distance(ti, "c1", "c2")
  d21 <- DropoutKiller:::.dk_tree_distance(ti, "c2", "c1")
  expect_true(is.finite(d12) && d12 >= 0 && d12 <= 1)
  expect_equal(d12, d21)
})

test_that("embedding-local weights decrease with biological projection distance", {
  z <- cbind(PC1 = c(0, 0.2, 1, 3, 6), PC2 = rep(0, 5))
  rownames(z) <- paste0("c", 1:5)
  membership <- setNames(rep(1L, 5), rownames(z))
  geom <- DropoutKiller:::.dk_build_local_geometry(
    z, membership, tree_weight = 0, local_k = 2, candidate_k = 4,
    min_effective_donors = 1, local_info_kappa = 1
  )
  w <- as.numeric(geom$W[1, ])
  expect_gt(w[2], w[4])
  expect_gt(w[4], w[5])
  expect_equal(w[1], 0)
  expect_equal(geom$W2, geom$W * geom$W)
  expect_equal(geom$total_weight, as.numeric(Matrix::rowSums(geom$W)))
})

test_that("hard biological strata remain absolute borrowing boundaries", {
  z <- cbind(PC1 = seq(0, 1, length.out = 8), PC2 = rep(0, 8))
  rownames(z) <- paste0("c", 1:8)
  membership <- setNames(rep(1L, 8), rownames(z))
  strata <- setNames(rep(c("A", "B"), each = 4), rownames(z))
  geom <- DropoutKiller:::.dk_build_local_geometry(
    z, membership, hard_stratum = strata, tree_weight = 0,
    local_k = 3, candidate_k = 7, min_effective_donors = 1,
    local_info_kappa = 1
  )
  expect_gt(sum(as.numeric(geom$W[1:4, 1:4])), 0)
  expect_gt(sum(as.numeric(geom$W[5:8, 5:8])), 0)
  expect_equal(sum(as.numeric(geom$W[1:4, 5:8])), 0)
  expect_equal(sum(as.numeric(geom$W[5:8, 1:4])), 0)
  expect_equal(geom$cell_stratum, rep(c("A", "B"), each = 4))
})

test_that("final membership is the recovery borrowing block even inside one hard stratum", {
  z <- cbind(PC1 = seq(0, 1, length.out = 8), PC2 = rep(0, 8))
  rownames(z) <- paste0("c", 1:8)
  membership <- setNames(rep(1:2, each = 4), rownames(z))
  strata <- setNames(rep("A", 8), rownames(z))
  geom <- DropoutKiller:::.dk_build_local_geometry(
    z, membership, hard_stratum = strata, tree_weight = 0,
    local_k = 3, candidate_k = 7, min_effective_donors = 1,
    local_info_kappa = 1
  )
  expect_gt(sum(as.numeric(geom$W[1:4, 1:4])), 0)
  expect_gt(sum(as.numeric(geom$W[5:8, 5:8])), 0)
  expect_equal(sum(as.numeric(geom$W[1:4, 5:8])), 0)
  expect_equal(sum(as.numeric(geom$W[5:8, 1:4])), 0)
  expect_equal(length(unique(geom$block_id)), 2L)
})

test_that("final memberships remain safe borrowing boundaries when no hard stratum is supplied", {
  set.seed(11)
  z <- cbind(PC1 = c(rnorm(6, -2, 0.1), rnorm(6, 2, 0.1)), PC2 = rnorm(12, 0, 0.05))
  rownames(z) <- paste0("c", 1:12)
  mf <- build_supercell_membership(z, gamma = 4, k_knn = 3,
                                   method = "walktrap", approximate = FALSE)
  expect_false(isTRUE(mf$settings$has_hard_stratum))
  geom <- DropoutKiller:::.dk_build_local_geometry(
    z, mf$membership, membership_fit = mf,
    tree_weight = 0.5, local_k = 3, candidate_k = 10,
    min_effective_donors = 1, local_info_kappa = 1
  )
  m <- mf$membership
  for (a in unique(m)) for (b in unique(m)) {
    block <- sum(as.numeric(geom$W[m == a, m == b, drop = FALSE]))
    if (a == b && sum(m == a) > 1) expect_gt(block, 0) else if (a != b) expect_equal(block, 0)
  }
})

test_that("batched local gene statistics equal scalar statistics", {
  set.seed(12)
  W <- matrix(runif(36), 6, 6); diag(W) <- 0
  W <- Matrix::Matrix(W, sparse = TRUE)
  geom <- list(W = W, W2 = W * W, total_weight = as.numeric(Matrix::rowSums(W)),
               Wt = Matrix::t(W), W2t = Matrix::t(W * W))
  x <- matrix(rexp(18), nrow = 3)
  x[cbind(c(1, 2, 3), c(1, 3, 5))] <- 0
  bat <- DropoutKiller:::.dk_local_gene_stats_batch(x, geom)
  for (g in seq_len(nrow(x))) {
    one <- DropoutKiller:::.dk_local_gene_stats(as.numeric(x[g, ]), geom)
    expect_equal(bat$mean[g, ], one$mean, tolerance = 1e-10)
    expect_equal(bat$variance[g, ], one$variance, tolerance = 1e-10)
    expect_equal(bat$effective_n[g, ], one$effective_n, tolerance = 1e-10)
    expect_equal(bat$prevalence[g, ], one$prevalence, tolerance = 1e-10)
  }
})

test_that("tree-local recovery uses a query-specific positive baseline", {
  cells <- paste0("c", 1:12)
  z <- cbind(PC1 = seq(0, 5.5, length.out = 12), PC2 = sin(seq_len(12)) / 20)
  rownames(z) <- cells
  x <- rbind(
    g1 = seq(1, 6.5, length.out = 12),
    g2 = seq(0.5, 3.25, length.out = 12),
    g3 = rep(c(1, 2), 6),
    g4 = seq(2, 3.1, length.out = 12),
    g5 = seq(0.3, 1.4, length.out = 12),
    g6 = seq(1.5, 2.6, length.out = 12)
  )
  colnames(x) <- cells
  truth <- x[1, 6]
  x[1, 6] <- 0
  mask <- Matrix::sparseMatrix(i = 1, j = 6, x = TRUE,
                               dims = dim(x), dimnames = dimnames(x))
  mf <- build_supercell_membership(z, gamma = 4, k_knn = 3,
                                   method = "walktrap", approximate = FALSE)
  rec <- recover_dropout_expression(
    x, mask, membership = mf, embedding = z,
    recovery_method = "tree_local_factor",
    min_target_observed = 100L,
    min_effective_donors = 1,
    tree_weight = 0.5, local_k = 3, candidate_k = 8,
    return_details = TRUE
  )
  expect_true(rec$events$changed[1])
  expect_equal(rec$events$recovery_method[1], "tree_local_mean")
  expect_equal(rec$events$recovered[1], rec$events$local_positive_mean[1])
  expect_true(is.finite(rec$events$effective_donors[1]) && rec$events$effective_donors[1] > 0)
  expect_true(is.finite(rec$events$prediction_sd[1]) && rec$events$prediction_sd[1] >= 0)
  local_stats <- DropoutKiller:::.dk_local_gene_stats(
    as.numeric(x[1, ]),
    list(W = rec$local_geometry$W, W2 = rec$local_geometry$W2,
         total_weight = rec$local_geometry$total_weight)
  )
  expect_equal(rec$events$effective_donors[1], local_stats$effective_n[6])
  expect_equal(rec$events$embedding_distance_weighted_mean[1],
               rec$local_geometry$embedding_distance_weighted_mean[6])
  expect_lt(abs(rec$events$recovered[1] - truth), max(abs(range(x[1, x[1, ] > 0]) - truth)))
})

test_that("tree-local mean fallback retains positive-expression outcome variance", {
  W <- matrix(1, 5, 5); diag(W) <- 0
  W <- Matrix::Matrix(W, sparse = TRUE)
  xg <- c(0, 1, 2, 4, 5)
  stats_g <- DropoutKiller:::.dk_local_gene_stats(xg, list(W = W))
  scores <- cbind(f1 = seq_len(5), f2 = seq_len(5)^2, f3 = c(1, 0, 1, 0, 1))
  fit <- DropoutKiller:::.dk_weighted_local_residual_target(
    xg, stats_g, scores, W, query = 1L,
    ridge = 1, min_target_observed = 2L,
    min_effective_donors = 1, local_info_kappa = 1
  )
  expect_equal(fit$method[1], "tree_local_mean")
  expect_true(is.finite(stats_g$variance[1]) && stats_g$variance[1] > 0)
  expect_gt(fit$prediction_sd[1]^2, stats_g$variance[1])
})

test_that("one gene-level residual model predicts multiple dropout queries", {
  set.seed(13)
  n <- 14
  W <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) W[i, j] <- exp(-abs(i - j) / 3)
  W <- Matrix::Matrix(W, sparse = TRUE)
  scores <- cbind(f1 = scale(seq_len(n))[, 1], f2 = scale(sin(seq_len(n) / 2))[, 1])
  xg <- 2 + 0.5 * scores[, 1] - 0.2 * scores[, 2]
  query <- c(4L, 10L); xg[query] <- 0
  stats_g <- DropoutKiller:::.dk_local_gene_stats(xg, list(W = W))
  fit <- DropoutKiller:::.dk_weighted_local_residual_target(
    xg, stats_g, scores, W, query = query,
    ridge = 1, min_target_observed = 5L,
    min_effective_donors = 1, local_info_kappa = 1
  )
  expect_length(fit$prediction, 2L)
  expect_true(all(is.finite(fit$prediction)))
  expect_true(all(fit$n_donors > 0))
  expect_true(all(fit$shrinkage >= 0 & fit$shrinkage <= 1))
})

test_that("tree-local ridge rejects invalid penalties", {
  W <- Matrix::Matrix(matrix(1, 5, 5) - diag(5), sparse = TRUE)
  xg <- c(0, 1, 2, 3, 4)
  stats_g <- DropoutKiller:::.dk_local_gene_stats(xg, list(W = W))
  scores <- matrix(seq_len(10), nrow = 5)
  expect_error(
    DropoutKiller:::.dk_weighted_local_residual_target(
      xg, stats_g, scores, W, query = 1L, ridge = -1,
      min_target_observed = 2L, min_effective_donors = 1
    ),
    "factor_ridge"
  )
})

test_that("tree-local engine preserves observed coordinates exactly", {
  set.seed(2)
  x <- matrix(rexp(8 * 20), nrow = 8)
  rownames(x) <- paste0("g", 1:8); colnames(x) <- paste0("c", 1:20)
  z <- cbind(PC1 = seq_len(20), PC2 = rnorm(20)); rownames(z) <- colnames(x)
  original <- x
  original[1, 10] <- 0; x <- original
  mask <- Matrix::sparseMatrix(i = 1, j = 10, x = TRUE,
                               dims = dim(x), dimnames = dimnames(x))
  mf <- build_supercell_membership(z, gamma = 5, k_knn = 3,
                                   method = "walktrap", approximate = FALSE)
  out <- recover_dropout_expression(
    x, mask, membership = mf, embedding = z,
    recovery_method = "tree_local_factor",
    min_target_observed = 3L, min_effective_donors = 1,
    candidate_k = 10L
  )
  keep <- matrix(TRUE, nrow(x), ncol(x)); keep[1, 10] <- FALSE
  expect_equal(as.numeric(out[keep]), as.numeric(x[keep]))
})