\name{build_supercell_membership}
\alias{build_supercell_membership}
\alias{dropout_membership}
\title{Build SuperCell-style cell memberships}
\description{Constructs biologically constrained memberships from a low-dimensional embedding using a Euclidean k-nearest-neighbor graph and graph community detection.}
\usage{
build_supercell_membership(embedding, group = NULL, split_by = NULL,
  gamma = 150, k_knn = 5L, method = c("walktrap", "louvain"),
  approximate = "auto", approx_n = 20000L, seed = 12345L,
  return_graph = FALSE)

dropout_membership(object = NULL, embedding = NULL, reduction = "pca",
  dims = 1:10, group = NULL, group_by = NULL, split_by = NULL,
  split_by_col = NULL, return_object = FALSE, ...)
}
\arguments{
\item{embedding}{Cell-by-dimension numeric matrix.}
\item{group}{Optional hard biological stratum, typically major cell type.}
\item{split_by}{Optional additional hard stratum such as condition or donor.}
\item{gamma}{Average graining level, default 150; target memberships are approximately \code{n/gamma}.}
\item{k_knn}{Number of nearest neighbors.}
\item{method}{Graph clustering method. Walktrap preserves the gamma-derived target.}
\item{approximate}{Logical or \code{"auto"}; use anchor/centroid approximation for large strata.}
\item{approx_n}{Approximate anchor count.}
\item{seed}{Random seed.}
\item{return_graph}{Whether to retain stratum graph objects.}
\item{object}{Optional Seurat object or embedding matrix for the convenience wrapper.}
\item{reduction}{Seurat reduction name.}
\item{dims}{Embedding dimensions to use.}
\item{group_by}{Seurat metadata column used as \code{group}.}
\item{split_by_col}{Seurat metadata column used as \code{split_by}.}
\item{return_object}{Return the full membership result rather than only the vector.}
\item{...}{Additional arguments passed to \code{build_supercell_membership}.}
}
\value{\code{build_supercell_membership} returns a \code{DropoutKillerMembership} object. \code{dropout_membership} returns the membership vector unless \code{return_object=TRUE}.}
