\name{supercell_lean_detect}
\alias{supercell_lean_detect}
\title{Detect recoverable zeros with SuperCell hierarchy and Lean membership features}
\description{Detect recoverable zeros with SuperCell hierarchy and Lean membership features.}
\usage{
supercell_lean_detect(x, embedding=NULL, group=NULL, model=NULL,
                                  geometry=NULL, normalize=TRUE, gamma=150,
                                  k_knn=5L, neighbor_k=30L, rank='auto', seed=1L)
}
\arguments{
\item{x}{Nonnegative gene-by-cell expression; raw counts if normalize is TRUE.}
\item{embedding}{RNA cell embedding, or NULL to derive it from the ALRA SVD.}
\item{group}{Broad cell classes. Fine subclusters are not inferred as boundaries.}
\item{model}{A model returned by fit_lean_detector.}
\item{geometry}{Optional precomputed internal geometry, including WNN affinity.}
\item{normalize}{Apply library-size normalization and log1p.}
\item{gamma}{Cells per metacell, with floor(n/gamma) Walktrap cut.}
\item{k_knn}{Graph neighbor count for RNA memberships.}
\item{neighbor_k}{Maximum support neighbors, borrowing only on shortage.}
\item{rank}{ALRA rank or auto.}
\item{seed}{Random seed for truncated SVD.}
}
\value{DropoutKillerDetection with called zero events and fitted model provenance.}
