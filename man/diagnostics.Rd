\name{validate_dropout_result}
\alias{validate_dropout_result}
\alias{membership_summary}
\title{DropoutKiller diagnostics}
\description{Checks selective-recovery invariants and summarizes membership size/purity.}
\usage{
validate_dropout_result(result, original, tolerance = 0)

membership_summary(membership, group = NULL)
}
\arguments{
\item{result}{A \code{DropoutKillerResult}.}
\item{original}{Original expression matrix.}
\item{tolerance}{Numerical comparison tolerance.}
\item{membership}{Membership vector.}
\item{group}{Optional biological group vector used to compute membership purity.}
}
\value{A named list of validation checks or a membership summary data frame.}
