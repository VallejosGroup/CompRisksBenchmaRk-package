\name{score_from_cifs}
\alias{score_from_cifs}
\title{Score CIF predictions against observed outcomes}
\usage{
score_from_cifs(
  out,
  test,
  times,
  causes,
  time_col = "time",
  status_col = "event",
  metrics = c("brier", "auc", "cindex_pec"),
  summary = c("ibs", "risks"),
  cens.method = "ipcw",
  cens.model = "km",
  cens.code = 0,
  se.fit = FALSE
)
}
\arguments{
  \item{out}{3-D array of CIF predictions with dimensions
    \code{[n, K, Tm]} where \code{n} is the number of test subjects,
    \code{K} the number of causes, and \code{Tm} the number of evaluation
    times.}
  \item{test}{Data frame of test observations containing at least the
    time and status columns.}
  \item{times}{Numeric vector of evaluation times (length \code{Tm}).}
  \item{causes}{Integer vector of cause codes (length \code{K}).}
  \item{time_col}{Name of the time column in \code{test}.}
  \item{status_col}{Name of the status/event column in \code{test}.}
  \item{metrics}{Character vector of metrics to compute.  Any subset of
    \code{"brier"}, \code{"auc"}, \code{"cidx_pec"}, \code{"cidx_survM"},
    \code{"calib_measures"}.}
  \item{summary}{Character vector; include \code{"ibs"} to also compute the
    integrated Brier score.}
  \item{cens.method}{Censoring weight method passed to
    \code{riskRegression::Score}.}
  \item{cens.model}{Censoring model; default \code{"km"}.}
  \item{cens.code}{Integer code for censored observations.}
  \item{se.fit}{Logical; compute standard errors (passed to
    \code{riskRegression::Score}).}
}
\description{
Iterates over causes and computes the requested performance metrics for each
cause separately using \code{riskRegression::Score} (Brier score / AUC),
\code{pec::cindex} (Uno-type C-index), and/or \code{\link{CalibrationPlot}}
(calibration measures).
}
\value{
A list with some or all of: \code{cause_bs}, \code{cause_ibs},
\code{cause_auc}, \code{cause_cindex_pec}, \code{cause_cindex_survM},
\code{cause_calib_measures}.  Each element is a list indexed by cause.
}
\seealso{
\code{\link{score_from_rmlt}}, \code{\link{CalibrationPlot}},
\code{\link{WeightedBrierScore}}
}
