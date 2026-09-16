\name{fit_lean_detector}
\alias{fit_lean_detector}
\title{Fit the balanced Lean hierarchy detector and a proxy-FPR threshold}
\description{Fit the balanced Lean hierarchy detector and a proxy-FPR threshold.}
\usage{
fit_lean_detector(features, truth, fpr=.01, method="Supercell_hierarchy_Lean_membership")
}
\arguments{
\item{features}{Data frame with z, local_hierarchy and H benchmark features.}
\item{truth}{Logical labels: artificially lost positives versus proxy negatives.}
\item{fpr}{Training proxy false-positive target, strictly between zero and one.}
\item{method}{RNA or WNN detector name; models are not interchangeable.}
}
\value{A reusable DropoutKillerLeanModel. Its score is not a biological posterior.}
