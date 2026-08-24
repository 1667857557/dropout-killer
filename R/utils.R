.dk_validate_expression <- function(x) {
  ok <- is.matrix(x) || inherits(x, "Matrix")
  if (!ok || length(dim(x)) != 2L) stop("x must be a matrix or Matrix object with genes in rows and cells in columns", call. = FALSE)
  if (nrow(x) < 1L || ncol(x) < 1L) stop("x must contain at least one gene and one cell", call. = FALSE)
  vals <- if (inherits(x, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(x)) methods::slot(x, "x") else rep.int(1, length(summary(x)$i))
  } else as.vector(x)
  if (!is.numeric(vals)) stop("x must contain numeric expression values", call. = FALSE)
  if (anyNA(vals) || (length(vals) && any(!is.finite(vals)))) stop("x contains non-finite values", call. = FALSE)
  if (length(vals) && any(vals < 0)) stop("x must be non-negative", call. = FALSE)
  x
}

.dk_names <- function(x) {
  genes <- rownames(x); cells <- colnames(x)
  if (is.null(genes)) genes <- paste0("gene_", seq_len(nrow(x)))
  if (is.null(cells)) cells <- paste0("cell_", seq_len(ncol(x)))
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) stop("gene names must be non-missing and unique", call. = FALSE)
  if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) stop("cell names must be non-missing and unique", call. = FALSE)
  list(genes = genes, cells = cells)
}

.dk_align_vector <- function(v, cells, name, allow_null = TRUE) {
  if (is.null(v)) {
    if (allow_null) return(NULL)
    stop(name, " is required", call. = FALSE)
  }
  if (!is.null(names(v))) {
    miss <- setdiff(cells, names(v))
    if (length(miss)) stop(name, " is missing ", length(miss), " cells", call. = FALSE)
    v <- v[cells]
  } else if (length(v) != length(cells)) stop(name, " must have one value per cell", call. = FALSE)
  unname(v)
}

.dk_align_embedding <- function(embedding, cells, dims = NULL) {
  if (is.null(embedding)) stop("embedding is required when membership is not supplied", call. = FALSE)
  z <- as.matrix(embedding)
  if (length(dim(z)) != 2L || ncol(z) < 1L) stop("embedding must be a cell-by-dimension matrix with at least one dimension", call. = FALSE)
  if (!is.numeric(z)) stop("embedding must be numeric", call. = FALSE)
  if (!is.null(rownames(z))) {
    miss <- setdiff(cells, rownames(z))
    if (length(miss)) stop("embedding is missing ", length(miss), " cells", call. = FALSE)
    z <- z[cells, , drop = FALSE]
  } else if (nrow(z) != length(cells)) stop("embedding must have one row per cell", call. = FALSE)
  if (!is.null(dims)) {
    dims <- as.integer(dims)
    dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(z)]
    if (!length(dims)) stop("no requested embedding dimensions are available", call. = FALSE)
    z <- z[, dims, drop = FALSE]
  }
  if (anyNA(z) || any(!is.finite(z))) stop("embedding contains non-finite values", call. = FALSE)
  rownames(z) <- cells
  z
}

.dk_align_membership <- function(membership, cells) {
  m <- .dk_align_vector(membership, cells, "membership", allow_null = FALSE)
  if (anyNA(m)) stop("membership contains NA values", call. = FALSE)
  as.integer(factor(m, levels = unique(m)))
}

.dk_stratum <- function(group, split_by, n) {
  clean <- function(v, default) {
    if (is.null(v)) return(rep(default, n))
    v <- as.character(v); v[is.na(v)] <- "<NA>"; v
  }
  g <- clean(group, "all")
  s <- clean(split_by, "all")
  interaction(g, s, drop = TRUE, lex.order = TRUE, sep = "||")
}

.dk_mask_events <- function(mask) {
  if (inherits(mask, "sparseMatrix")) {
    sm <- summary(mask)
    keep <- if ("x" %in% names(sm)) as.logical(sm$x) else rep(TRUE, nrow(sm))
    return(data.frame(i = sm$i[keep], j = sm$j[keep]))
  }
  if (!is.matrix(mask)) mask <- as.matrix(mask)
  idx <- which(mask, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(i = integer(), j = integer()))
  data.frame(i = idx[, 1L], j = idx[, 2L])
}

.dk_sparse_numeric <- function(i, j, value, nr, nc, dimnames = NULL) {
  if (!length(i)) return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(nr, nc), dimnames = dimnames))
  Matrix::sparseMatrix(i = i, j = j, x = value, dims = c(nr, nc), dimnames = dimnames)
}

.dk_sparse_logical <- function(i, j, nr, nc, dimnames = NULL) {
  if (!length(i)) return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(nr, nc), dimnames = dimnames) > 0)
  Matrix::sparseMatrix(i = i, j = j, x = TRUE, dims = c(nr, nc), dimnames = dimnames)
}

.dk_row_means <- function(x) {
  if (inherits(x, "Matrix")) Matrix::rowMeans(x) else rowMeans(x)
}

.dk_row_sums <- function(x) {
  if (inherits(x, "Matrix")) Matrix::rowSums(x) else rowSums(x)
}
