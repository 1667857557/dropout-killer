.dk_alra_library_log <- function(x, scale_factor = 1e4) {
  x <- .dk_validate_expression(x)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L || !is.finite(scale_factor) || scale_factor <= 0)
    stop("normalization_scale_factor must be a finite value > 0", call. = FALSE)
  lib <- if (inherits(x, "Matrix")) Matrix::colSums(x) else colSums(x)
  bad <- !is.finite(lib) | lib <= 0
  if (any(bad))
    stop("ALRA library+log normalization requires a positive library size for every cell; remove zero-library cells before DropoutKiller so expression, embedding, and membership remain aligned", call. = FALSE)
  sf <- as.numeric(scale_factor / lib)
  if (inherits(x, "Matrix")) {
    out <- x %*% Matrix::Diagonal(n = length(sf), x = sf)
    if (inherits(out, "sparseMatrix") && "x" %in% methods::slotNames(out)) {
      methods::slot(out, "x") <- log1p(methods::slot(out, "x"))
    } else {
      out <- log1p(out)
    }
  } else {
    out <- sweep(x, 2L, lib, FUN = "/")
    out <- log1p(out * scale_factor)
  }
  dimnames(out) <- dimnames(x)
  attr(out, "dropoutkiller_normalization") <- list(
    method = "ALRA_library_size_log1p",
    scale_factor = as.numeric(scale_factor)
  )
  out
}
