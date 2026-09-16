\name{build_supercell_membership}
\alias{build_supercell_membership}
\title{Build SuperCell-style cell memberships}
\description{Constructs biologically constrained memberships from a low-dimensional embedding using a Euclidean k-nearest-neighbor graph. Walktrap retains the full hierarchy so downstream recovery can use graded tree proximity.}
\usage{
build_supercell_membership(embedding, group = NULL, split_by = NULL,
  gamma = 150, k_knn = 5L, method = c("walktrap", "louvain"),
  approximate = "auto", approx_n = 20000L, seed = 12345L,
  return_graph = FALSE)
}
\arguments{
\item{embedding}{Cell-by-dimension numeric matrix.}
\item{group}{Optional hard biological stratum, typically a broad cell lineage.}
\item{split_by}{Optional additional hard stratum for direct geometry construction. The high-level Lean detector itself uses broad group boundaries only.}
\item{gamma}{Average graining level, default 150; target memberships are approximately \code{n/gamma}.}
\item{k_knn}{Number of nearest neighbors used to build each stratum graph.}
\item{method}{Graph clustering method. Walktrap provides a retained hierarchy; Louvain has no hierarchy.}
\item{approximate}{Logical or \code{"auto"}; use anchor/centroid approximation for large strata.}
\item{approx_n}{Approximate anchor count.}
\item{seed}{Random seed.}
\item{return_graph}{Whether to additionally retain stratum graph objects.}
}
\value{A \code{DropoutKillerMembership} object containing \code{membership}, \code{membership_table}, per-stratum summaries, named \code{cell_stratum}, retained Walktrap \code{hierarchies}, cached \code{tree_indices}, and optional graphs.}
