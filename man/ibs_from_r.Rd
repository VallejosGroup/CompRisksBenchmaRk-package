\name{ibs_from_r}
\alias{ibs_from_r}
\title{Integrated Brier score wrapper for Python interop}
\usage{
ibs_from_r(
  out,
  e_va,
  t_va,
  times,
  causes,
  time_col = "time",
  status_col = "event",
  cens.method = "ipcw",
  cens.model = "km",
  cens.code = 0,
  force_sequential = TRUE,
  quiet = TRUE,
  log_file = NULL,
  debug = FALSE,
  debug_n = 5
)
}
\arguments{
  \item{out}{3-D array \code{[n, K, Tm]} of CIF predictions.}
  \item{e_va}{Integer vector of event codes for the validation set.}
  \item{t_va}{Numeric vector of observed times for the validation set.}
  \item{times}{Numeric vector of evaluation times (length \code{Tm}).}
  \item{causes}{Integer vector of cause codes (length \code{K}).}
  \item{time_col}{Name to use for the time column.}
  \item{status_col}{Name to use for the status column.}
  \item{cens.method}{Censoring method passed to \code{riskRegression::Score}.}
  \item{cens.model}{Censoring model.}
  \item{cens.code}{Integer censoring code.}
  \item{force_sequential}{Logical; force single-threaded execution to avoid
    parallelism conflicts when called from Python.}
  \item{quiet}{Logical; suppress output.}
  \item{log_file}{Optional path to a log file.}
  \item{debug}{Logical; return extra debugging information.}
  \item{debug_n}{Number of rows to include in debug output.}
}
\description{
Computes the IBS per cause using \code{riskRegression::Score}.  Designed to
be called from Python (via \code{rpy2}) where the CIF array is produced by
a Python model.
}
\value{
A numeric matrix \code{[K, Tm]} of IBS values (one row per cause, one column
per evaluation time), or a list when \code{debug = TRUE}.
}
