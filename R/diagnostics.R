#' Summarize membership size and purity
#'
#' @export
membership_summary <- function(membership, group = NULL) {
  m <- as.integer(factor(membership, levels = unique(membership)))
  tab <- as.data.frame(table(m), stringsAsFactors = FALSE); names(tab) <- c("membership", "n_cells")
  tab$membership <- as.integer(as.character(tab$membership))
  if (!is.null(group)) {
    if (length(group) != length(m)) stop("group and membership lengths differ", call. = FALSE)
    g <- as.character(group); g[is.na(g)] <- "<NA>"
    purity <- vapply(sort(unique(m)), function(k) max(table(g[m == k])) / sum(m == k), numeric(1))
    tab$group_purity <- purity[match(tab$membership, sort(unique(m)))]
  }
  tab
}

#' Validate zero-selective recovery invariants
#'
#' @export
validate_dropout_result <- function(result, original, tolerance = 0) {
  if (!inherits(result, "DropoutKillerResult")) stop("result must be a DropoutKillerResult", call. = FALSE)
  original <- .dk_validate_expression(original)
  same_dim <- identical(dim(result$expression), dim(original))
  if (!same_dim) return(list(valid = FALSE, same_dimensions = FALSE, observed_nonzero_preserved = FALSE, only_masked_zeros_changed = FALSE, nonnegative = FALSE))
  if (inherits(original, "sparseMatrix")) {
    sm <- Matrix::summary(original); ov <- if (nrow(sm)) as.numeric(original[cbind(sm$i, sm$j)]) else numeric()
    nv <- if (nrow(sm)) as.numeric(result$expression[cbind(sm$i, sm$j)]) else numeric()
  } else {
    idx <- which(original != 0, arr.ind = TRUE); ov <- if (nrow(idx)) original[idx] else numeric(); nv <- if (nrow(idx)) result$expression[idx] else numeric()
  }
  preserved <- length(ov) == length(nv) && all(abs(ov - nv) <= tolerance)
  ev <- .dk_mask_events(result$mask)
  masked_zero <- if (nrow(ev)) all(as.vector(original[cbind(ev$i, ev$j)] == 0)) else TRUE
  if (inherits(result$expression, "sparseMatrix") || inherits(original, "sparseMatrix")) {
    delta <- result$expression - original; ds <- Matrix::summary(delta)
    changed_keys <- if (nrow(ds)) paste(ds$i[abs(ds$x) > tolerance], ds$j[abs(ds$x) > tolerance], sep = ":") else character()
  } else {
    di <- which(abs(result$expression - original) > tolerance, arr.ind = TRUE)
    changed_keys <- if (nrow(di)) paste(di[, 1L], di[, 2L], sep = ":") else character()
  }
  mask_keys <- if (nrow(ev)) paste(ev$i, ev$j, sep = ":") else character()
  only_mask <- all(changed_keys %in% mask_keys) && masked_zero
  vals <- if (inherits(result$expression, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(result$expression)) methods::slot(result$expression, "x") else rep.int(1, length(Matrix::summary(result$expression)$i))
  } else as.vector(result$expression)
  nonnegative <- !anyNA(vals) && all(vals >= -tolerance)
  checks <- list(same_dimensions = same_dim, observed_nonzero_preserved = preserved,
                 only_masked_zeros_changed = only_mask, nonnegative = nonnegative)
  c(list(valid = all(unlist(checks))), checks)
}
