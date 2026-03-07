\name{CindexCR}
\alias{CindexCR}
\title{Concordance index for competing risks (Wolbers)}
\usage{
CindexCR(time, status, predicted, Cause_int = 1)
}
\arguments{
  \item{time}{Numeric vector of observed times (must be positive).}
  \item{status}{Integer vector of event codes (0, 1, 2).}
  \item{predicted}{Numeric vector of predicted CIFs for \code{Cause_int}.}
  \item{Cause_int}{Integer cause of interest (default 1).}
}
\description{
Computes the concordance index for competing risks following the notation of
Wolbers et al.  Evaluation is performed at the maximum observed time.
}
\value{
A list with \code{cindex} (numeric) and \code{time_ev} (the evaluation
time, i.e. \code{max(time) + 1}).
}
\references{
Wolbers M et al. (2009). Prognostic models with competing risks: methods and
application to coronary risk prediction. \emph{Epidemiology}, 20(4):555-561.
}
