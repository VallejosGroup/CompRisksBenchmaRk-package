\name{IntegratedBrierScore}
\alias{IntegratedBrierScore}
\title{Integrated Brier score for competing risks}
\usage{
IntegratedBrierScore(
  prediction_matrix,
  taus,
  time,
  status,
  cause,
  cens.code,
  cmprsk
)
}
\arguments{
  \item{prediction_matrix}{Numeric matrix \code{[n, length(taus)]} of
    predicted CIFs.}
  \item{taus}{Numeric vector of evaluation times.}
  \item{time}{Numeric vector of observed times.}
  \item{status}{Integer vector of event codes.}
  \item{cause}{Integer code for the event of interest.}
  \item{cens.code}{Integer code for censored observations.}
  \item{cmprsk}{Logical; passed to \code{\link{WeightedBrierScore}}.}
}
\description{
Integrates the IPCW Brier score over time using the trapezoidal rule,
normalising by the time interval \code{[min(taus), max(taus)]}.
}
\value{
A list with \code{integrated.brier.score}, \code{average.brier.score},
\code{weighted.brier.score}, and \code{taus}.
}
\seealso{
\code{\link{WeightedBrierScore}}
}
