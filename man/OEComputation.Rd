\name{OEComputation}
\alias{OEComputation}
\title{Observed-to-expected ratio with confidence interval}
\usage{
OEComputation(data, pred_mat, cause, times, horizon)
}
\arguments{
  \item{data}{Data frame with \code{time} and \code{status} columns.}
  \item{pred_mat}{Numeric matrix \code{[n, Tm]} of predicted CIFs.}
  \item{cause}{Integer cause code.}
  \item{times}{Numeric vector of evaluation times (length \code{Tm}).}
  \item{horizon}{Scalar evaluation time at which to compute the O/E ratio.}
}
\description{
Computes the observed-to-expected (O/E) ratio at a given horizon time using
the Aalen-Johansen estimator for the observed cumulative incidence and the
mean predicted CIF as the expected value.  A confidence interval is derived
via the delta method.
}
\value{
A formatted string (invisibly) and side-effect: a named numeric vector
\code{c(OE, lower, upper)}.
}
\seealso{
\code{\link{CalibrationPlot}}
}
