test_that("membership never crosses hard groups", {
  set.seed(1)
  z <- rbind(matrix(rnorm(40, 0, 0.1), 20, 2), matrix(rnorm(40, 5, 0.1), 20, 2))
  rownames(z) <- paste0("c", seq_len(nrow(z)))
  g <- rep(c("A", "B"), each = 20)
  fit <- build_supercell_membership(z, group = g, gamma = 10, k_knn = 3, seed = 7)
  pur <- membership_summary(fit$membership, g)
  expect_true(all(pur$group_purity == 1))
  expect_equal(length(fit$membership), nrow(z))
})

test_that("membership is reproducible", {
  set.seed(2); z <- matrix(rnorm(120), 40, 3)
  a <- build_supercell_membership(z, gamma = 10, k_knn = 4, seed = 9)$membership
  b <- build_supercell_membership(z, gamma = 10, k_knn = 4, seed = 9)$membership
  expect_equal(a, b)
})
