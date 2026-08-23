#' Build a target-by-source gene network from an edge list
#'
#' The returned orientation follows the prediction contract A[target, source].
#' For a GRN edge TF -> target, use the TF column as source and the regulated
#' gene column as target. Set directed = FALSE for symmetric PPI-style edges.
#'
#' @export
gene_network_from_edges <- function(edges, genes = NULL, source = "source", target = "target",
                                    weight = NULL, directed = TRUE) {
  if (!is.data.frame(edges)) stop("edges must be a data.frame", call. = FALSE)
  if (!source %in% names(edges) || !target %in% names(edges)) stop("source/target columns were not found", call. = FALSE)
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
  if (!length(i)) return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(length(genes), length(genes)), dimnames = list(genes, genes)))
  Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(length(genes), length(genes)),
                       dimnames = list(genes, genes), use.last.ij = FALSE)
}

#' Prepare a gene-network adjacency matrix
#'
#' @export
prepare_gene_network <- function(adjacency, genes, remove_self = TRUE, normalize = TRUE) {
  if (is.null(adjacency)) return(NULL)
  if (!(is.matrix(adjacency) || inherits(adjacency, "Matrix"))) stop("adjacency must be a matrix or Matrix object", call. = FALSE)
  genes <- as.character(genes); G <- length(genes)
  if (!G || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) stop("genes must be non-empty, non-missing, and unique", call. = FALSE)
  A <- Matrix::Matrix(adjacency, sparse = TRUE)
  sm0 <- summary(A)
  if ("x" %in% names(sm0) && nrow(sm0) && any(!is.finite(as.numeric(sm0$x)))) stop("adjacency contains non-finite values", call. = FALSE)
  rn <- rownames(A); cn <- colnames(A)
  if (is.null(rn) || is.null(cn)) {
    if (!identical(dim(A), c(G, G))) stop("unnamed adjacency must have dimensions length(genes) x length(genes)", call. = FALSE)
    dimnames(A) <- list(genes, genes)
  } else {
    if (anyDuplicated(rn) || anyDuplicated(cn)) stop("adjacency row and column names must be unique", call. = FALSE)
    sm <- summary(A)
    ri <- match(rn[sm$i], genes); cj <- match(cn[sm$j], genes)
    keep <- !is.na(ri) & !is.na(cj) & is.finite(sm$x)
    A <- Matrix::sparseMatrix(i = ri[keep], j = cj[keep], x = sm$x[keep], dims = c(G, G), dimnames = list(genes, genes))
  }
  if (remove_self) diag(A) <- 0
  A <- Matrix::drop0(A)
  if (normalize) {
    rs <- Matrix::rowSums(abs(A)); inv <- ifelse(rs > 0, 1 / rs, 0)
    A <- Matrix::Diagonal(x = inv) %*% A
  }
  Matrix::drop0(A)
}

#' Combine PPI, pathway, GRN, or other gene priors
#'
#' @export
combine_gene_prior <- function(networks, genes, weights = NULL, remove_self = TRUE) {
  if (is.null(networks)) return(NULL)
  if (!is.list(networks) || inherits(networks, "Matrix") || is.matrix(networks)) networks <- list(networks)
  keep <- !vapply(networks, is.null, logical(1))
  if (!any(keep)) return(NULL)
  if (!is.null(weights)) {
    if (length(weights) == length(networks)) weights <- weights[keep]
    else if (length(weights) != sum(keep)) stop("weights must match all supplied networks or all non-NULL networks", call. = FALSE)
  }
  networks <- networks[keep]
  if (is.null(weights)) weights <- rep(1 / length(networks), length(networks))
  if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
    stop("weights must be non-negative finite values matching networks", call. = FALSE)
  weights <- weights / sum(weights)
  mats <- lapply(networks, prepare_gene_network, genes = genes, remove_self = remove_self, normalize = TRUE)
  A <- mats[[1L]] * weights[1L]
  if (length(mats) > 1L) for (i in 2:length(mats)) A <- A + mats[[i]] * weights[i]
  rs <- Matrix::rowSums(abs(A)); A <- Matrix::Diagonal(x = ifelse(rs > 0, 1 / rs, 0)) %*% A
  Matrix::drop0(A)
}

.dk_gene_predict_events <- function(x, adjacency, membership, events, standardize = TRUE, min_sd = 1e-8) {
  n_ev <- nrow(events); pred <- rep(NA_real_, n_ev); support <- integer(n_ev)
  if (!n_ev || is.null(adjacency)) return(list(prediction = pred, support = support))
  lev <- unique(events$membership)
  for (m in lev) {
    cells <- which(membership == m); if (!length(cells)) next
    ev_m <- which(events$membership == m)
    targets <- sort(unique(events$i[ev_m]))
    W <- adjacency[targets, , drop = FALSE]
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

#' Gene-network prediction at masked events
#'
#' @export
gene_prior_prediction <- function(x, adjacency, membership, mask, standardize = TRUE,
                                  min_sd = 1e-8, return_events = FALSE) {
  x <- .dk_validate_expression(x); nm <- .dk_names(x)
  membership <- .dk_align_membership(membership, nm$cells)
  A <- prepare_gene_network(adjacency, nm$genes)
  if (length(dim(mask)) != 2L || !identical(as.integer(dim(mask)), as.integer(dim(x)))) stop("mask and x dimensions differ", call. = FALSE)
  events <- .dk_mask_events(mask)
  if (nrow(events)) events$membership <- membership[events$j] else events$membership <- integer()
  fit <- .dk_gene_predict_events(x, A, membership, events, standardize, min_sd)
  if (return_events) {
    events$prediction <- fit$prediction; events$network_support <- fit$support
    return(events)
  }
  ok <- is.finite(fit$prediction)
  .dk_sparse_numeric(events$i[ok], events$j[ok], fit$prediction[ok], nrow(x), ncol(x), list(nm$genes, nm$cells))
}
