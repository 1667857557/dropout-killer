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
})

test_that("hard biological strata remain absolute borrowing boundaries", {
  z <- cbind(PC1 = seq(0, 1, length.out = 8), PC2 = rep(0, 8))
  rownames(z) <- paste0("c", 1:8)
  membership <- setNames(rep(1:2, each = 4), rownames(z))
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

test_that("tree-local recovery uses a query-specific positive baseline", {
  cells <- paste0("c", 1:12)
  genes <- paste0("g", 1:6)
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
  expect_lt(abs(rec$events$recovered[1] - truth), max(abs(range(x[1, x[1, ] > 0]) - truth)))
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
