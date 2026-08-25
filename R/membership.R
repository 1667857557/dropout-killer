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

.dk_index_walktrap_tree <- function(hierarchy, leaf_names) {
  if (is.null(hierarchy) || !igraph::is_hierarchical(hierarchy)) return(NULL)
  merge <- tryCatch(igraph::merges(hierarchy), error = function(e) NULL)
  n <- length(leaf_names)
  if (is.null(merge) || !n || !nrow(merge)) return(NULL)
  merge <- as.matrix(merge)
  max_node <- n + nrow(merge)
  parent <- integer(max_node)
  subtree_size <- rep.int(1L, max_node)
  for (r in seq_len(nrow(merge))) {
    node <- n + r
    child <- as.integer(merge[r, ])
    child <- child[is.finite(child) & child >= 1L & child < node]
    if (!length(child)) next
    parent[child] <- node
    subtree_size[node] <- sum(subtree_size[child])
  }
  denom <- log(max(2, n))
  height <- numeric(max_node)
  internal <- seq.int(n + 1L, max_node)
  height[internal] <- pmin(1, log(pmax(2, subtree_size[internal])) / denom)
  ancestors <- vector("list", n)
  for (i in seq_len(n)) {
    path <- i; cur <- i
    guard <- 0L
    while (parent[cur] > 0L && guard <= max_node) {
      cur <- parent[cur]
      path <- c(path, cur)
      guard <- guard + 1L
    }
    ancestors[[i]] <- path
  }
  names(ancestors) <- leaf_names
  list(parent = parent, subtree_size = subtree_size, height = height,
       ancestors = ancestors, leaf_names = leaf_names, n_leaves = n)
}

.dk_tree_distance <- function(tree_index, query_names, donor_names) {
  if (length(query_names) != length(donor_names)) stop("query_names and donor_names must have equal length", call. = FALSE)
  out <- rep(NA_real_, length(query_names))
  if (is.null(tree_index) || !length(query_names)) return(out)
  a <- tree_index$ancestors
  for (k in seq_along(query_names)) {
    qi <- a[[query_names[k]]]; dj <- a[[donor_names[k]]]
    if (is.null(qi) || is.null(dj)) next
    common <- qi[qi %in% dj]
    if (!length(common)) { out[k] <- 1; next }
    lca <- common[1L]
    out[k] <- tree_index$height[lca]
  }
  out
}

.dk_cluster_graph <- function(z, target_n, k_knn, method = "walktrap") {
  n <- nrow(z)
  if (n <= 1L || target_n >= n) {
    return(list(membership = seq_len(n), graph = .dk_knn_graph(z, k_knn), hierarchy = NULL))
  }
  g <- .dk_knn_graph(z, k_knn)
  ncomp <- igraph::components(g)$no
  target_n <- max(as.integer(target_n), ncomp)
  if (target_n >= n) return(list(membership = seq_len(n), graph = g, hierarchy = NULL))
  hierarchy <- NULL
  if (method == "walktrap") {
    hierarchy <- igraph::cluster_walktrap(g, merges = TRUE)
    mem <- igraph::cut_at(hierarchy, no = target_n)
  } else if (method == "louvain") {
    warning("method='louvain' ignores gamma-derived target membership count and has no hierarchy", call. = FALSE)
    mem <- igraph::membership(igraph::cluster_louvain(g))
  } else stop("method must be 'walktrap' or 'louvain'", call. = FALSE)
  list(membership = as.integer(factor(mem, levels = unique(mem))), graph = g, hierarchy = hierarchy)
}

.dk_cluster_stratum <- function(z, gamma, k_knn, method, approximate, approx_n, seed) {
  n <- nrow(z)
  target_n <- max(1L, min(n, as.integer(round(n / gamma))))
  use_approx <- isTRUE(approximate) || (identical(approximate, "auto") && n > 50000L)
  if (!use_approx || n <= approx_n || target_n >= n) {
    out <- .dk_cluster_graph(z, target_n, k_knn, method)
    return(list(membership = out$membership, graph = out$graph, hierarchy = out$hierarchy,
                hierarchy_rows = seq_len(n), approximate = FALSE, target_n = target_n, anchor_n = n))
  }
  anchor_n <- min(n, max(as.integer(approx_n), 3L * target_n))
  if (anchor_n >= n) {
    out <- .dk_cluster_graph(z, target_n, k_knn, method)
    return(list(membership = out$membership, graph = out$graph, hierarchy = out$hierarchy,
                hierarchy_rows = seq_len(n), approximate = FALSE, target_n = target_n, anchor_n = n))
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
  list(membership = as.integer(factor(mem, levels = unique(mem))), graph = fit$graph,
       hierarchy = fit$hierarchy, hierarchy_rows = anchor, approximate = TRUE,
       target_n = target_n, anchor_n = anchor_n)
}

#' Build SuperCell-style memberships from a low-dimensional embedding
#'
#' Cells are first partitioned into hard biological strata (for example major
#' cell type, and optionally condition/donor), then each stratum is represented
#' by a Euclidean kNN graph. Walktrap hierarchical clustering is cut to roughly
#' n/gamma memberships. The full walktrap hierarchy is retained so recovery can
#' use continuous tree proximity instead of treating the final cut as an
#' equal-weight bag of cells.
#'
#' @export
build_supercell_membership <- function(embedding, group = NULL, split_by = NULL, gamma = 150,
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
  has_hard_stratum <- !is.null(group) || !is.null(split_by)
  group <- .dk_align_vector(group, cells, "group")
  split_by <- .dk_align_vector(split_by, cells, "split_by")
  strata <- .dk_stratum(group, split_by, n)
  membership <- integer(n); info <- list(); graphs <- list(); hierarchies <- list(); tree_indices <- list(); offset <- 0L
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
    hierarchies[[lev[s]]] <- fit$hierarchy
    hcells <- cells[idx[fit$hierarchy_rows]]
    tree_indices[[lev[s]]] <- .dk_index_walktrap_tree(fit$hierarchy, hcells)
  }
  membership <- as.integer(factor(membership, levels = unique(membership)))
  names(membership) <- cells
  cell_stratum <- as.character(strata); names(cell_stratum) <- cells
  tab <- as.data.frame(table(membership), stringsAsFactors = FALSE)
  names(tab) <- c("membership", "n_cells"); tab$membership <- as.integer(as.character(tab$membership))
  out <- list(membership = membership, membership_table = tab, strata = do.call(rbind, info),
              cell_stratum = cell_stratum, hierarchies = hierarchies, tree_indices = tree_indices,
              settings = list(gamma = gamma, k_knn = as.integer(k_knn), method = method,
                              approximate = approximate, approx_n = as.integer(approx_n), seed = as.integer(seed),
                              has_hard_stratum = has_hard_stratum))
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
      if (!split_by_col %in% colnames(meta)) stop("split_by_col column not found in Seurat metadata", call. = FALSE)
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
  if (length(x$hierarchies)) cat(" hierarchical strata:", sum(vapply(x$hierarchies, function(h) !is.null(h), logical(1))), "\n")
  invisible(x)
}
