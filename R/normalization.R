.dk_prepare_expression_layers <- function(counts, scale_factor = 1e4) {
  counts <- .dk_validate_expression(counts)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0)
    stop("scale_factor must be a finite value > 0", call. = FALSE)
  lib <- if (inherits(counts, "Matrix")) Matrix::colSums(counts) else colSums(counts)
  if (any(!is.finite(lib) | lib <= 0))
    stop("raw count layers require a positive library size for every cell", call. = FALSE)
  sf <- as.numeric(scale_factor / lib)
  if (inherits(counts, "Matrix")) {
    log_expression <- counts %*% Matrix::Diagonal(n = length(sf), x = sf)
    if (inherits(log_expression, "sparseMatrix") &&
        "x" %in% methods::slotNames(log_expression)) {
      methods::slot(log_expression, "x") <- log1p(methods::slot(log_expression, "x"))
    } else log_expression <- log1p(log_expression)
  } else {
    log_expression <- log1p(sweep(counts, 2L, lib, "/") * scale_factor)
  }
  dimnames(log_expression) <- dimnames(counts)
  attr(log_expression, "dropoutkiller_normalization") <- list(
    method = "ALRA_library_size_log1p", scale_factor = as.numeric(scale_factor)
  )
  list(
    counts = counts,
    library_size = as.numeric(lib),
    mean_library_size = mean(lib),
    size_factor = as.numeric(lib / mean(lib)),
    log_expression = log_expression,
    scale_factor = as.numeric(scale_factor)
  )
}

.dk_alra_library_log <- function(x, scale_factor = 1e4) {
  x <- .dk_validate_expression(x)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L || !is.finite(scale_factor) || scale_factor <= 0)
    stop("normalization_scale_factor must be a finite value > 0", call. = FALSE)
  .dk_prepare_expression_layers(x, scale_factor)$log_expression
}
