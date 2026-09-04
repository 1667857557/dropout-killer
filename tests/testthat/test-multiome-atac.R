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

test_that("small WNN strata stay on the supplied weighted graph when approximation is requested", {
  cells <- paste0("c", seq_len(12))
  a <- Matrix::Matrix(0, 12, 12, sparse = TRUE)
  for (i in 1:5) { a[i, i + 1] <- 0.4 + i / 10; a[i + 1, i] <- a[i, i + 1] }
  for (i in 7:11) { a[i, i + 1] <- 0.4 + (i - 6) / 10; a[i + 1, i] <- a[i, i + 1] }
  rownames(a) <- colnames(a) <- cells
  group <- rep(c("A", "B"), each = 6); names(group) <- cells
  fit <- build_supercell_membership(graph = a, group = group, gamma = 3,
    method = "walktrap", approximate = TRUE, approx_n = 10,
    graph_dims = 3, return_graph = TRUE)
  expect_false(any(fit$strata$approximate))
  expect_true(all(fit$strata$geometry == "weighted_graph"))
  expect_equal(sort(igraph::E(fit$graphs[["A"]])$weight), sort(a[1:6, 1:6][a[1:6, 1:6] > 0]) / 2,
    tolerance = 1e-12)
})

test_that("pattern sparse WNN adjacency is accepted", {
  a <- Matrix::Matrix(0, 6, 6, sparse = TRUE)
  for (i in 1:5) { a[i, i + 1] <- 1; a[i + 1, i] <- 1 }
  rownames(a) <- colnames(a) <- paste0("c", 1:6)
  pattern <- methods::as(a > 0, "nMatrix")
  expect_false("x" %in% methods::slotNames(pattern))
  g <- DropoutKiller:::.dk_as_weighted_graph(pattern)
  expect_true(igraph::is_igraph(g))
  expect_equal(igraph::vcount(g), 6)
})

test_that("WNN spectral embedding uses leading algebraic eigenmodes", {
  n <- 8L
  a <- Matrix::Matrix(0, n, n, sparse = TRUE)
  for (i in seq_len(n - 1L)) { a[i, i + 1L] <- 1; a[i + 1L, i] <- 1 }
  rownames(a) <- colnames(a) <- paste0("c", seq_len(n))
  z <- DropoutKiller:::.dk_graph_spectral_embedding(a, dims = 2)
  deg <- as.numeric(Matrix::rowSums(a))
  d <- 1 / sqrt(pmax(deg, .Machine$double.eps))
  s <- Matrix::Diagonal(x = d) %*% a %*% Matrix::Diagonal(x = d)
  ee <- eigen(as.matrix(s), symmetric = TRUE)
  use <- order(ee$values, decreasing = TRUE)[1:2]
  target <- ee$vectors[, use, drop = FALSE]
  expect_lt(max(abs(z %*% t(z) - target %*% t(target))), 1e-5)
})

test_that("ATAC binarization respects stored and dense zero values", {
  sparse_zero <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(0, 2), dims = c(2, 2))
  dense_zero <- Matrix::Matrix(matrix(c(0, 0, 0, 2), 2, 2), sparse = FALSE)
  expected <- matrix(c(0, 0, 0, 1), 2, 2)
  expect_equal(as.matrix(DropoutKiller:::.dk_atac_binary_matrix(sparse_zero)), expected)
  expect_equal(as.matrix(DropoutKiller:::.dk_atac_binary_matrix(dense_zero)), expected)
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

test_that("Seurat wrappers reject recovered assay names that overwrite inputs", {
  skip_if_not_installed("Seurat")
  counts <- matrix(c(1, 0, 0, 1), 2, 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2")))
  object <- Seurat::CreateSeuratObject(counts = counts)
  expect_error(dropout_killer_seurat(object, assay = "RNA", modality = "rna", new_assay = "RNA"),
    "must differ from the input assay")

  object[["ATAC"]] <- Seurat::CreateAssayObject(counts = counts)
  expect_error(dropout_killer_seurat(object, assay = "RNA", modality = "multiome",
    atac_assay = "ATAC", rna_new_assay = "RNA", atac_new_assay = "ATAC_DK"),
    "must differ from both raw RNA and ATAC")
  expect_error(dropout_killer_seurat(object, assay = "RNA", modality = "multiome",
    atac_assay = "ATAC", rna_new_assay = "DK", atac_new_assay = "DK"),
    "must be different")
})
