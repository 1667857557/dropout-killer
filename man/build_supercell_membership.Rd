\name{build_supercell_membership}
\alias{build_supercell_membership}
\alias{dropout_membership}
\title{Build SuperCell-style cell memberships}
\description{Constructs biologically constrained memberships from either a low-dimensional embedding or a weighted multimodal graph. For weighted WNN graphs (for example a Seurat \code{wsnn} graph), the graph itself drives the Walktrap cut and a normalized-adjacency spectral embedding is retained for downstream hierarchy/geometry-aware recovery. The full Walktrap hierarchy is retained whenever available.}
\usage{
build_supercell_membership(embedding = NULL, group = NULL, split_by = NULL,
  gamma = 150, k_knn = 5L, method = c("walktrap", "louvain"),
  approximate = "auto", approx_n = 20000L, seed = 12345L,
  return_graph = FALSE, graph = NULL, graph_dims = 30L)

dropout_membership(object = NULL, embedding = NULL, reduction = "pca",
  dims = 1:10, group = NULL, group_by = NULL, split_by = NULL,
  split_by_col = NULL, return_object = FALSE, graph = NULL,
  graph_name = NULL, graph_dims = 30L, ...)
}
\arguments{
\item{embedding}{Optional cell-by-dimension numeric matrix. Required when \code{graph=NULL}.}
\item{group}{Optional hard biological stratum, typically major cell type.}
\item{split_by}{Optional additional hard stratum such as condition or donor.}
\item{gamma}{Average graining level, default 150; target memberships are approximately \code{n/gamma}.}
\item{k_knn}{Number of nearest neighbors used when an embedding-derived graph is required.}
\item{method}{Graph clustering method. Walktrap retains a hierarchy; Louvain has no hierarchy.}
\item{approximate}{Logical or \code{"auto"}; use anchor/centroid approximation for large strata. Large WNN strata use the retained WNN spectral embedding for approximation.}
\item{approx_n}{Approximate anchor count.}
\item{seed}{Random seed.}
\item{return_graph}{Whether to additionally retain stratum graph objects.}
\item{graph}{Optional weighted graph supplied as an \code{igraph} object or square matrix-like object.}
\item{graph_dims}{Number of normalized-adjacency spectral dimensions retained for recovery geometry.}
\item{object}{Optional Seurat object or embedding matrix for the convenience wrapper.}
\item{reduction}{Seurat reduction name when \code{graph_name} is not supplied.}
\item{dims}{Embedding dimensions to use.}
\item{group_by}{Seurat metadata column used as \code{group}.}
\item{split_by_col}{Seurat metadata column used as \code{split_by}.}
\item{return_object}{Return the full membership result rather than only the vector.}
\item{graph_name}{Optional graph stored in \code{object@graphs}, for example \code{"wsnn"}.}
\item{...}{Additional arguments passed to \code{build_supercell_membership}.}
}
\value{\code{build_supercell_membership} returns a \code{DropoutKillerMembership} object containing membership assignments, retained hierarchy/tree indices, and a recovery embedding. \code{dropout_membership} returns only the membership vector unless \code{return_object=TRUE}.}
