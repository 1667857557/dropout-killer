test_that("simplex projection returns non-negative unit-sum weights", {
  w <- DropoutKiller:::.dk_project_simplex(c(-2, 0.2, 1.7, 0.4))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-12)
})

test_that("barycentric geometry treats membership as a soft prior", {
  cells <- paste0("c", 1:8)
  z <- cbind(PC1 = c(0, 0.1, 0.2, 0.3, 0.31, 0.4, 0.5, 0.6), PC2 = rep(0, 8))
  rownames(z) <- cells
  membership <- setNames(rep(1:2, each = 4), cells)
  geom <- DropoutKiller:::.dk_build_barycentric_geometry(
    z, membership, query_cells = 4L, hard_stratum = setNames(rep("A", 8), cells),
    tree_weight = 0, local_k = 3, candidate_k = 7,
    membership_penalty = 0.5, barycentric_lambda = 1,
    barycentric_iter = 30L, barycentric_tol = 1e-8,
    min_effective_donors = 1
  )
  w <- as.numeric(geom$W[1, ])
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_true(all(w >= 0))
  expect_gt(sum(w[5:8]), 0)
  expect_true(is.finite(geom$state_error[1]))
  expect_true(is.finite(geom$prior_error[1]))
})

test_that("barycentric recovery can borrow target expression across nearby memberships", {
  cells <- paste0("c", 1:8)
  z <- cbind(PC1 = c(0, 0.1, 0.2, 0.3, 0.31, 0.4, 0.5, 0.6), PC2 = rep(0, 8))
  rownames(z) <- cells
  x <- rbind(g1 = c(0, 0, 0, 0, 5, 6, 7, 8),
             g2 = seq(1, 1.7, length.out = 8),
             g3 = seq(2, 2.7, length.out = 8))
  colnames(x) <- cells
  membership <- setNames(rep(1:2, each = 4), cells)
  strata <- setNames(rep("A", 8), cells)
  mask <- Matrix::sparseMatrix(i = 1, j = 4, x = TRUE, dims = dim(x), dimnames = dimnames(x))
  out <- recover_dropout_expression(
    x, mask, membership = membership, embedding = z, hard_stratum = strata,
    recovery_method = "barycentric", tree_weight = 0,
    local_k = 3, candidate_k = 7, membership_penalty = 0.25,
    barycentric_lambda = 1, min_effective_donors = 1,
    return_details = TRUE
  )
  expect_true(out$events$changed[1])
  expect_equal(out$events$recovery_method[1], "barycentric_mean")
  expect_gt(out$events$recovered[1], 0)
  expect_true(is.finite(out$events$prediction_sd[1]))
  keep <- matrix(TRUE, nrow(x), ncol(x)); keep[1, 4] <- FALSE
  expect_equal(as.numeric(out$expression[keep]), as.numeric(x[keep]))
})

test_that("hard strata remain absolute boundaries for barycentric recovery", {
  cells <- paste0("c", 1:8)
  z <- cbind(PC1 = c(0, 0.1, 0.2, 0.3, 0.31, 0.4, 0.5, 0.6), PC2 = rep(0, 8))
  rownames(z) <- cells
  x <- rbind(g1 = c(0, 0, 0, 0, 5, 6, 7, 8),
             g2 = seq(1, 1.7, length.out = 8))
  colnames(x) <- cells
  membership <- setNames(rep(1:2, each = 4), cells)
  strata <- setNames(rep(c("A", "B"), each = 4), cells)
  mask <- Matrix::sparseMatrix(i = 1, j = 4, x = TRUE, dims = dim(x), dimnames = dimnames(x))
  out <- recover_dropout_expression(
    x, mask, membership = membership, embedding = z, hard_stratum = strata,
    recovery_method = "barycentric", tree_weight = 0,
    local_k = 3, candidate_k = 7, membership_penalty = 0,
    barycentric_lambda = 1, min_effective_donors = 1,
    return_details = TRUE
  )
  expect_false(out$events$changed[1])
  expect_equal(out$expression[1, 4], 0)
  expect_equal(sum(as.numeric(out$local_geometry$W[1, 5:8])), 0)
})

test_that("barycentric recovery reports finite uncertainty when support is adequate", {
  cells <- paste0("c", 1:12)
  z <- cbind(PC1 = seq(0, 1.1, length.out = 12), PC2 = sin(seq_len(12)) / 20)
  rownames(z) <- cells
  x <- rbind(g1 = seq(1, 6.5, length.out = 12),
             g2 = seq(0.5, 3.25, length.out = 12),
             g3 = rep(c(1, 2), 6))
  colnames(x) <- cells
  x[1, 6] <- 0
  membership <- setNames(rep(1:3, each = 4), cells)
  mask <- Matrix::sparseMatrix(i = 1, j = 6, x = TRUE, dims = dim(x), dimnames = dimnames(x))
  out <- recover_dropout_expression(
    x, mask, membership = membership, embedding = z,
    recovery_method = "barycentric", tree_weight = 0,
    local_k = 4, candidate_k = 10, membership_penalty = 0.25,
    barycentric_lambda = 2, min_effective_donors = 3,
    return_details = TRUE
  )
  expect_true(out$events$changed[1])
  expect_true(out$uncertainty_available)
  expect_true(is.finite(out$events$prediction_sd[1]) && out$events$prediction_sd[1] >= 0)
  expect_true(out$events$effective_donors[1] >= 3)
})
