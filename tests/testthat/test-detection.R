test_that("selection can never overwrite observed non-zero values", {
  x <- matrix(c(5, 0, 4, 3, 0, 2), nrow = 2)
  score <- matrix(0.99, nrow(x), ncol(x))
  mask <- select_dropout_mask(score, threshold = 0.95, x = x)
  idx <- which(as.matrix(mask), arr.ind = TRUE)
  expect_true(all(x[idx] == 0))
})

test_that("local ALRA score is sparse and zero-only", {
  set.seed(3)
  x <- matrix(rexp(120), 12, 10); x[x < 0.5] <- 0
  colnames(x) <- paste0("c", 1:10); rownames(x) <- paste0("g", 1:12)
  s <- local_alra_score(x, rep(1, 10), rank = 2, quantile_prob = 0.1, min_cells = 4)
  expect_true(inherits(s, "sparseMatrix"))
  sm <- summary(s)
  if (nrow(sm)) expect_true(all(x[cbind(sm$i, sm$j)] == 0))
})
