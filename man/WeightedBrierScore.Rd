\name{WeightedBrierScore}
\alias{WeightedBrierScore}
\title{IPCW-weighted Brier score for competing risks}
\usage{
WeightedBrierScore(
  predictions,
  tau,
  time,
  status,
  cause,
  cens.code,
  cmprsk = FALSE
)
}
\arguments{
  \item{predictions}{Numeric matrix \code{[n, length(tau)]} of predicted
    CIFs, or a vector when \code{length(tau) == 1}.}
  \item{tau}{Numeric vector of evaluation times.}
  \item{time}{Numeric vector of observed times.}
  \item{status}{Integer vector of event codes.}
  \item{cause}{Integer code for the event of interest.}
  \item{cens.code}{Integer code for censored observations (typically 0).}
  \item{cmprsk}{Logical; if \code{TRUE} applies the competing-risks
    weighting scheme; otherwise uses standard survival weighting.}
}
\description{
Computes the inverse probability of censoring weighted (IPCW) Brier score at
one or more evaluation times.  Censoring probabilities are estimated via a
marginal Kaplan-Meier model (treating all non-censored events as events).
}
\value{
A list with elements \code{weighted.brier.score} (numeric vector, one value
per element of \code{tau}), \code{tau}, \code{n}, and \code{n.risk}.
}
\seealso{
\code{\link{IntegratedBrierScore}}
}
