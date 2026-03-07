\name{compute_rmlt_from_cif}
\alias{compute_rmlt_from_cif}
\title{Compute restricted mean lifetime from CIF}
\usage{
compute_rmlt_from_cif(cif, times, tau)
}
\arguments{
  \item{cif}{3-D array \code{[n, K, Tm]} of cumulative incidence functions.}
  \item{times}{Numeric vector of evaluation times (length \code{Tm}).}
  \item{tau}{Scalar upper truncation time for the integral.}
}
\description{
Integrates each subject's CIF curve up to \code{tau} using the trapezoidal
rule to obtain the restricted mean lifetime lost (RMLT) due to each cause.
}
\value{
A matrix \code{[n, K]} of RMLT values.
}
\seealso{
\code{\link{score_from_rmlt}}
}
