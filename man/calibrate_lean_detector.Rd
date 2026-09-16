\name{calibrate_lean_detector}
\alias{calibrate_lean_detector}
\title{Calibrate Lean hierarchy detection using artificial UMI thinning}
\description{Calibrate Lean hierarchy detection using artificial UMI thinning.}
\usage{
calibrate_lean_detector(counts, group, q=.5, seeds=10001:10003,
                                    fpr=.01, negative_max=.005, negative_other_min=.2,
                                    gamma=150, k_knn=5L, neighbor_k=30L, rank='auto',
                                    geometry_builder=NULL)
}
\arguments{
\item{counts}{Raw integer RNA counts (genes by cells).}
\item{group}{Broad cell classes used to define proxy negatives and hard boundaries.}
\item{q}{UMI retention probability, not the fraction of matrix entries masked.}
\item{seeds}{Independent thinning seeds. All resulting events are used.}
\item{fpr}{Training proxy-FPR target.}
\item{negative_max}{Maximum prevalence in the negative's own broad class.}
\item{negative_other_min}{Minimum prevalence in another broad class.}
\item{gamma}{Cells per SuperCell membership.}
\item{k_knn}{RNA graph neighbor count.}
\item{neighbor_k}{Maximum hierarchy support neighbors.}
\item{rank}{ALRA rank or auto.}
\item{geometry_builder}{Optional function(counts) returning freshly rebuilt WNN geometry.}
}
\value{DropoutKillerLeanModel; use separate masks/data for evaluation.}
