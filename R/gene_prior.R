#' Build a PPI or pathway prior from an edge list
#'
#' The returned orientation follows A[target, source]. For PPI, symmetric edges
#' are usually appropriate; pathway edges may be directed when the source data
#' encode directionality. GRN/TF-target priors are intentionally unsupported.
#'
#' @export
gene_prior_from_edges <- function(edges, genes = NULL, source = "source", target = "target",
                                  weight = NULL, prior_type = c("ppi", "pathway"),
                                  directed = NULL) {
  prior_type <- match.arg(prior_type)
  if (!is.data.frame(edges)) stop("edges must be a data.frame", call. = FALSE)
  if (!source %in% names(edges) || !target %in% names(edges)) stop("source/target columns were not found", call. = FALSE)
  if (is.null(directed)) directed <- identical(prior_type, "pathway")
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) stop("directed must be TRUE or FALSE", call. = FALSE)
  src <- as.character(edges[[source]]); tgt <- as.character(edges[[target]])
  w <- if (is.null(weight)) rep(1, nrow(edges)) else {
    if (!weight %in% names(edges)) stop("weight column was not found", call. = FALSE)
    as.numeric(edges[[weight]])
  }
  keep <- !is.na(src) & nzchar(src) & !is.na(tgt) & nzchar(tgt) & is.finite(w)
  src <- src[keep]; tgt <- tgt[keep]; w <- w[keep]
  if (is.null(genes)) genes <- unique(c(tgt, src))
  genes <- as.character(genes)
  if (!length(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) stop("genes must be non-empty, non-missing, and unique", call. = FALSE)
  i <- match(tgt, genes); j <- match(src, genes); keep <- !is.na(i) & !is.na(j)
  i <- i[keep]; j <- j[keep]; w <- w[keep]
  if (!isTRUE(directed)) { oi <- i; oj <- j; ow <- w; i <- c(oi, oj); j <- c(oj, oi); w <- c(ow, ow) }
  if (!length(i)) {
    A <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(length(genes), length(genes)),
                              dimnames = list(genes, genes))
  } else {
    A <- Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(length(genes), length(genes)),
                              dimnames = list(genes, genes), use.last.ij = FALSE)
  }
  attr(A, "prior_type") <- prior_type
  A
}

#' Prepare an aligned PPI or pathway prior
#'
#' @export
prepare_gene_prior <- function(prior, genes, remove_self = TRUE, normalize = TRUE) {
  if (is.null(prior)) return(NULL)
  if (!(is.matrix(prior) || inherits(prior, "Matrix"))) stop("prior must be a matrix or Matrix object", call. = FALSE)
  genes <- as.character(genes); G <- length(genes)
  if (!G || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) stop("genes must be non-empty, non-missing, and unique", call. = FALSE)
  A <- Matrix::Matrix(prior, sparse = TRUE)
  sm0 <- summary(A)
  if ("x" %in% names(sm0) && nrow(sm0) && any(!is.finite(as.numeric(sm0$x)))) stop("prior contains non-finite values", call. = FALSE)
  rn <- rownames(A); cn <- colnames(A)
  if (is.null(rn) || is.null(cn)) {
    if (!identical(dim(A), c(G, G))) stop("unnamed prior must have dimensions length(genes) x length(genes)", call. = FALSE)
    dimnames(A) <- list(genes, genes)
  } else {
    if (anyDuplicated(rn) || anyDuplicated(cn)) stop("prior row and column names must be unique", call. = FALSE)
    sm <- summary(A)
    ri <- match(rn[sm$i], genes); cj <- match(cn[sm$j], genes)
    keep <- !is.na(ri) & !is.na(cj) & is.finite(sm$x)
    A <- Matrix::sparseMatrix(i = ri[keep], j = cj[keep], x = sm$x[keep],
                              dims = c(G, G), dimnames = list(genes, genes))
  }
  if (remove_self) diag(A) <- 0
  A <- Matrix::drop0(A)
  if (normalize) {
    rs <- Matrix::rowSums(abs(A)); inv <- ifelse(rs > 0, 1 / rs, 0)
    A <- Matrix::Diagonal(x = inv) %*% A
  }
  Matrix::drop0(A)
}

#' Combine PPI and pathway priors
#'
#' Only PPI and pathway inputs are accepted by the public workflow. If both are
#' supplied, each prior is aligned and row-normalized before weighted fusion.
#'
#' @export
combine_gene_prior <- function(genes, ppi = NULL, pathway = NULL, weights = NULL, remove_self = TRUE) {
  present <- c(ppi = !is.null(ppi), pathway = !is.null(pathway))
  if (!any(present)) return(NULL)
  priors <- list(ppi = ppi, pathway = pathway)[present]
  if (is.null(weights)) {
    weights <- rep(1 / length(priors), length(priors)); names(weights) <- names(priors)
  } else {
    if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) stop("weights must be non-negative finite numeric values", call. = FALSE)
    if (!is.null(names(weights))) {
      if (any(!names(weights) %in% c("ppi", "pathway"))) stop("named weights may contain only ppi and pathway", call. = FALSE)
      if (!all(names(priors) %in% names(weights))) stop("named weights must contain all supplied priors", call. = FALSE)
      weights <- weights[names(priors)]
    } else if (length(weights) == 2L) {
      names(weights) <- c("ppi", "pathway"); weights <- weights[names(priors)]
    } else if (length(weights) != length(priors)) {
      stop("weights must have length 2 (ppi, pathway) or match the supplied priors", call. = FALSE)
    }
    if (sum(weights) <= 0) stop("weights must sum to a positive value", call. = FALSE)
    weights <- weights / sum(weights)
  }
  mats <- lapply(priors, prepare_gene_prior, genes = genes, remove_self = remove_self, normalize = TRUE)
  A <- mats[[1L]] * weights[1L]
  if (length(mats) > 1L) for (i in 2:length(mats)) A <- A + mats[[i]] * weights[i]
  rs <- Matrix::rowSums(abs(A))
  A <- Matrix::Diagonal(x = ifelse(rs > 0, 1 / rs, 0)) %*% A
  Matrix::drop0(A)
}

.dk_prior_predict_events <- function(x, prior, membership, events, standardize = TRUE, min_sd = 1e-8) {
  n_ev <- nrow(events); pred <- rep(NA_real_, n_ev); support <- integer(n_ev)
  if (!n_ev || is.null(prior)) return(list(prediction = pred, support = support))
  lev <- unique(events$membership)
  for (m in lev) {
    cells <- which(membership == m); if (!length(cells)) next
    ev_m <- which(events$membership == m)
    targets <- sort(unique(events$i[ev_m]))
    W <- prior[targets, , drop = FALSE]
    support_target <- Matrix::rowSums(abs(W) > 0)
    target_map <- match(events$i[ev_m], targets)
    support[ev_m] <- support_target[target_map]
    active <- support_target > 0
    if (!any(active)) next
    xk <- x[, cells, drop = FALSE]
    cand_cells <- sort(unique(events$j[ev_m])); xq <- x[, cand_cells, drop = FALSE]
    if (standardize) {
      mu <- .dk_row_means(xk)
      ex2 <- .dk_row_means(xk * xk)
      if (length(cells) > 1L) vv <- pmax(0, (ex2 - mu * mu) * length(cells) / (length(cells) - 1L)) else vv <- rep(0, length(mu))
      sdv <- sqrt(vv); inv_sd <- ifelse(is.finite(sdv) & sdv > min_sd, 1 / sdv, 0)
      B <- W %*% Matrix::Diagonal(x = inv_sd)
      baseline <- as.vector(W %*% (mu * inv_sd))
      zp <- as.matrix(B %*% xq)
      zp <- sweep(zp, 1L, baseline, "-")
      tp <- sweep(zp, 1L, sdv[targets], "*")
      tp <- sweep(tp, 1L, mu[targets], "+")
      valid_target <- active & is.finite(sdv[targets]) & sdv[targets] > min_sd
    } else {
      tp <- as.matrix(W %*% xq); valid_target <- active
    }
    tp[!is.finite(tp)] <- NA_real_; tp[tp < 0] <- 0
    col_map <- match(events$j[ev_m], cand_cells)
    vals <- tp[cbind(target_map, col_map)]
    ok <- valid_target[target_map] & is.finite(vals)
    pred[ev_m[ok]] <- vals[ok]
  }
  list(prediction = pred, support = support)
}

#' PPI/pathway prior prediction at masked events
#'
#' @export
gene_prior_prediction <- function(x, prior, membership, mask, standardize = TRUE,
                                  min_sd = 1e-8, return_events = FALSE) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  A <- prepare_gene_prior(prior, nm$genes)
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x)))) stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) events$membership <- membership[events$j] else events$membership <- integer()
  fit <- .dk_prior_predict_events(x, A, membership, events, standardize, min_sd)
  if (return_events) {
    events$prediction <- fit$prediction; events$prior_support <- fit$support
    return(events)
  }
  ok <- is.finite(fit$prediction)
  .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok], nrow(x), ncol(x), list(nm$genes, nm$cells))
}
