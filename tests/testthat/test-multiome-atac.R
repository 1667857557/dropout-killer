test_that("weighted graph memberships retain a WNN-derived geometry", {
  cells <- paste0("c", seq_len(12))
  a <- Matrix::Matrix(0, 12, 12, sparse = TRUE)
  for (i in 1:5) { a[i, i + 1] <- 1; a[i + 1, i] <- 1 }
  for (i in 7:11) { a[i, i + 1] <- 1; a[i + 1, i] <- 1 }
  rownames(a) <- colnames(a) <- cells
  group <- rep(c("A", "B"), each = 6); names(group) <- cells
  fit <- build_supercell_membership(graph = a, group = group, gamma = 3,
    method = "walktrap", approximate = FALSE, graph_dims = 3)
  expect_s3_class(fit, "DropoutKillerMembership")
  expect_equal(rownames(fit$embedding), cells)
  expect_equal(fit$settings$geometry_source, "weighted_graph")
  pur <- membership_summary(fit$membership, group)
  expect_true(all(pur$group_purity == 1))
})

test_that("ATAC detector and recovery only alter selected zero coordinates", {
  cells <- paste0("c", seq_len(10)); peaks <- paste0("p", seq_len(20))
  x <- matrix(0, 20, 10, dimnames = list(peaks, cells))
  x[1:12, 1:8] <- 1; x[1:2, 9:10] <- 1
  x[3, 1:8] <- c(1, 2, 1, 1, 2, 1, 1, 2)
  x <- Matrix::Matrix(x, sparse = TRUE)
  z <- cbind(WNN1 = seq(0, 0.09, length.out = 10), WNN2 = rep(0, 10))
  rownames(z) <- cells; membership <- setNames(rep(1L, 10), cells)
  fit <- DropoutKiller:::.dk_atac_dropout_killer(counts = x, embedding = z,
    membership = membership, bfdr = 0.20, pi_min = 0.60, pre_pi = 0.40,
    score_min = 0.10, min_observed_donors = 2, min_effective_donors = 2,
    candidate_k = 9, local_k = 5, phi_prior = 1)
  expect_true(nrow(fit$events) > 0)
  expect_true(all(x[cbind(fit$events$i, fit$events$j)] == 0))
  observed <- which(x != 0, arr.ind = TRUE)
  expect_equal(as.numeric(fit$expression[observed]), as.numeric(x[observed]))
  expect_true(all(fit$events$dropout_probability >= 0 & fit$events$dropout_probability <= 1))
  expect_true(all(fit$events$recovered[fit$events$changed] > 0))
})

test_that("ATAC capture correction follows the PIC product model", {
  cells <- paste0("c", seq_len(8)); peaks <- paste0("p", seq_len(12))
  x <- matrix(0, 12, 8, dimnames = list(peaks, cells)); x[, 1:6] <- 1; x[1:2, 7:8] <- 1
  x <- Matrix::Matrix(x, sparse = TRUE)
  cap <- DropoutKiller:::.dk_atac_capture_probability(x, max_iter = 50, tol = 1e-7)$capture
  expect_true(all(is.finite(cap))); expect_true(all(cap > 0 & cap < 1))
  expect_gt(mean(cap[1:6]), mean(cap[7:8]))
})
