.dk_knn_graph <- function(z, k_knn) {
  n <- nrow(z)
  if (n <= 1L) return(igraph::make_empty_graph(n = n, directed = FALSE))
  k_use <- min(as.integer(k_knn), n - 1L)
  if (k_use < 1L) stop("k_knn must be >= 1", call. = FALSE)
  nn <- RANN::nn2(data = z, query = z, k = min(n, k_use + 1L))$nn.idx
  from <- integer(); to <- integer()
  for (i in seq_len(n)) {
    v <- nn[i, ]; v <- v[v != i]
    if (length(v) > k_use) v <- v[seq_len(k_use)]
    if (length(v)) { from <- c(from, rep.int(i, length(v))); to <- c(to, v) }
  }
  if (!length(from)) return(igraph::make_empty_graph(n = n, directed = FALSE))
  g <- igraph::graph_from_edgelist(cbind(from, to), directed = FALSE)
  if (igraph::vcount(g) < n) g <- igraph::add_vertices(g, n - igraph::vcount(g))
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

.dk_cluster_graph <- function(z, target_n, k_knn, method = "walktrap") {
  n <- nrow(z)
  if (n <= 1L || target_n >= n) return(list(membership = seq_len(n), graph = .dk_knn_graph(z, k_knn)))
  g <- .dk_knn_graph(z, k_knn)
  ncomp <- igraph::components(g)$no
  target_n <- max(as.integer(target_n), ncomp)
  if (target_n >= n) return(list(membership = seq_len(n), graph = g))
  if (method == "walktrap") {
    cl <- igraph::cluster_walktrap(g)
    mem <- igraph::cut_at(cl, no = target_n)
  } else if (method == "louvain") {
    warning("method='louvain' ignores gamma-derived target membership count", call. = FALSE)
    mem <- igraph::membership(igraph::cluster_louvain(g))
  } else stop("method must be 'walktrap' or 'louvain'", call. = FALSE)
  list(membership = as.integer(factor(mem, levels = unique(mem))), graph = g)
}

.dk_cluster_stratum <- function(z, gamma, k_knn, method, approximate, approx_n, seed) {
  n <- nrow(z)
  target_n <- max(1L, min(n, as.integer(round(n / gamma))))
  use_approx <- isTRUE(approximate) || (identical(approximate, "auto") && n > 50000L)
  if (!use_approx || n <= approx_n || target_n >= n) {
    out <- .dk_cluster_graph(z, target_n, k_knn, method)
    return(list(membership = out$membership, graph = out$graph, approximate = FALSE, target_n = target_n, anchor_n = n))
  }
  anchor_n <- min(n, max(as.integer(approx_n), 3L * target_n))
  if (anchor_n >= n) {
    out <- .dk_cluster_graph(z, target_n, k_knn, method)
    return(list(membership = out$membership, graph = out$graph, approximate = FALSE, target_n = target_n, anchor_n = n))
  }
  set.seed(seed)
  anchor <- sort(sample.int(n, anchor_n, replace = FALSE))
  rest <- setdiff(seq_len(n), anchor)
  fit <- .dk_cluster_graph(z[anchor, , drop = FALSE], target_n, k_knn, method)
  am <- fit$membership
  lev <- sort(unique(am))
  centroids <- vapply(lev, function(k) colMeans(z[anchor[am == k], , drop = FALSE]), numeric(ncol(z)))
  centroids <- t(centroids)
  assigned <- RANN::nn2(data = centroids, query = z[rest, , drop = FALSE], k = 1L)$nn.idx[, 1L]
  mem <- integer(n); mem[anchor] <- am; mem[rest] <- lev[assigned]
  list(membership = as.integer(factor(mem, levels = unique(mem))), graph = fit$graph, approximate = TRUE, target_n = target_n, anchor_n = anchor_n)
}

#' Build SuperCell-style memberships from a low-dimensional embedding
#'
#' Cells are first partitioned into hard biological strata (for example major
#' cell type, and optionally condition/donor), then each stratum is represented
#' by a Euclidean kNN graph. Walktrap hierarchical clustering is cut to roughly
#' n/gamma memberships. This is intentionally stricter than post-hoc splitting:
#' graph construction itself cannot borrow edges across supplied strata.
#'
#' @export
build_supercell_membership <- function(embedding, group = NULL, split_by = NULL, gamma = 20,
                                       k_knn = 5L, method = c("walktrap", "louvain"),
                                       approximate = "auto", approx_n = 20000L,
                                       seed = 12345L, return_graph = FALSE) {
  z <- as.matrix(embedding)
  if (length(dim(z)) != 2L || nrow(z) < 1L || ncol(z) < 1L) stop("embedding must be a non-empty cell-by-dimension matrix", call. = FALSE)
  if (!is.numeric(z)) stop("embedding must be numeric", call. = FALSE)
  if (anyNA(z) || any(!is.finite(z))) stop("embedding contains non-finite values", call. = FALSE)
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 2) stop("gamma must be a finite number >= 2", call. = FALSE)
  if (!is.numeric(k_knn) || length(k_knn) != 1L || !is.finite(k_knn) || k_knn < 1) stop("k_knn must be >= 1", call. = FALSE)
  if (!is.numeric(approx_n) || length(approx_n) != 1L || !is.finite(approx_n) || approx_n < 2) stop("approx_n must be >= 2", call. = FALSE)
  method <- match.arg(method)
  if (!((is.logical(approximate) && length(approximate) == 1L && !is.na(approximate)) || identical(approximate, "auto"))) stop("approximate must be TRUE, FALSE, or 'auto'", call. = FALSE)
  n <- nrow(z); cells <- rownames(z); if (is.null(cells)) cells <- paste0("cell_", seq_len(n))
  if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) stop("embedding row names must be non-missing and unique", call. = FALSE)
  group <- .dk_align_vector(group, cells, "group")
  split_by <- .dk_align_vector(split_by, cells, "split_by")
  strata <- .dk_stratum(group, split_by, n)
  membership <- integer(n); info <- list(); graphs <- list(); offset <- 0L
  lev <- levels(strata)
  for (s in seq_along(lev)) {
    idx <- which(strata == lev[s])
    fit <- .dk_cluster_stratum(z[idx, , drop = FALSE], gamma, k_knn, method, approximate, approx_n, seed + s - 1L)
    local <- fit$membership
    membership[idx] <- local + offset
    offset <- offset + max(local)
    info[[s]] <- data.frame(stratum = lev[s], n_cells = length(idx), target_memberships = fit$target_n,
                            observed_memberships = length(unique(local)), approximate = fit$approximate,
                            anchor_n = fit$anchor_n, stringsAsFactors = FALSE)
    if (return_graph) graphs[[lev[s]]] <- fit$graph
  }
  names(membership) <- cells
  tab <- as.data.frame(table(membership), stringsAsFactors = FALSE)
  names(tab) <- c("membership", "n_cells"); tab$membership <- as.integer(as.character(tab$membership))
  out <- list(membership = membership, membership_table = tab, strata = do.call(rbind, info),
              settings = list(gamma = gamma, k_knn = as.integer(k_knn), method = method,
                              approximate = approximate, approx_n = as.integer(approx_n), seed = as.integer(seed)))
  if (return_graph) out$graphs <- graphs
  class(out) <- "DropoutKillerMembership"
  out
}

#' SuperCell-style membership convenience wrapper
#'
#' @export
dropout_membership <- function(object = NULL, embedding = NULL, reduction = "pca", dims = 1:10,
                               group = NULL, group_by = NULL, split_by = NULL, split_by_col = NULL,
                               return_object = FALSE, ...) {
  if (!is.null(object) && inherits(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required for a Seurat object", call. = FALSE)
    embedding <- Seurat::Embeddings(object, reduction = reduction)
    if (!is.null(dims)) {
      dims <- as.integer(dims); dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(embedding)]
      if (!length(dims)) stop("no requested embedding dimensions are available", call. = FALSE)
      embedding <- embedding[, dims, drop = FALSE]
    }
    meta <- object[[]]
    if (!is.null(group_by)) {
      if (!group_by %in% colnames(meta)) stop("group_by column not found in Seurat metadata", call. = FALSE)
      group <- meta[[group_by]]; names(group) <- rownames(meta)
    }
    if (!is.null(split_by_col)) {
      if (!split_by_col %in% colnames(meta)) stop("split_by_col not found in Seurat metadata", call. = FALSE)
      split_by <- meta[[split_by_col]]; names(split_by) <- rownames(meta)
    }
  } else if (is.null(embedding) && !is.null(object)) embedding <- object
  if (is.null(embedding)) stop("provide embedding or a Seurat object", call. = FALSE)
  fit <- build_supercell_membership(embedding = embedding, group = group, split_by = split_by, ...)
  if (return_object) fit else fit$membership
}

#' @export
print.DropoutKillerMembership <- function(x, ...) {
  cat("DropoutKiller membership\n")
  cat(" cells:", length(x$membership), "\n")
  cat(" memberships:", length(unique(x$membership)), "\n")
  cat(" median size:", stats::median(x$membership_table$n_cells), "\n")
  invisible(x)
}
