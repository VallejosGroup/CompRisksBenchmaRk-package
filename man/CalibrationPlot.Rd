\name{CalibrationPlot}
\alias{CalibrationPlot}
\title{Calibration plot for competing risks predictions}
\usage{
CalibrationPlot(
  model_name = NULL,
  predictions,
  data,
  time,
  status,
  tau,
  cause = 1,
  cens.code = 0,
  loess_smoothing = TRUE,
  predictions.type = c("CIF", "survival"),
  bandwidth = NULL,
  graph = TRUE
)
}
\arguments{
  \item{model_name}{Optional character label used in plot titles.}
  \item{predictions}{Numeric matrix \code{[n, length(tau)]} or vector (when
    \code{length(tau) == 1}) of predicted CIFs (or survival probabilities).}
  \item{data}{Data frame used to compute pseudovalues.}
  \item{time}{Numeric vector of observed times.}
  \item{status}{Integer vector of event codes.}
  \item{tau}{Numeric vector of evaluation times.}
  \item{cause}{Integer cause code.}
  \item{cens.code}{Integer censoring code.}
  \item{loess_smoothing}{Logical; if \code{TRUE} (default) smoothing is done
    via LOESS; otherwise \code{prodlim::meanNeighbors} (nearest-neighbour) is
    used.}
  \item{predictions.type}{One of \code{"CIF"} or \code{"survival"};
    survival predictions are converted to CIF via \code{1 - predictions}.}
  \item{bandwidth}{Optional bandwidth for nearest-neighbour smoothing
    (ignored when \code{loess_smoothing = TRUE}).}
  \item{graph}{Logical; if \code{TRUE} ggplot2 calibration graphs are
    returned.}
}
\description{
Produces calibration plots based on jackknife pseudovalues (Aalen-Johansen
estimator via \code{prodlim}) smoothed with LOESS or nearest-neighbour
kernel smoothing.  Numerical calibration measures (ICI, E50, E90, Emax, RSB)
and observed/expected (O/E) ratio are returned alongside the plots.

Inspired by Van Geloven et al. (2022).
}
\value{
A list with:
\describe{
  \item{graphs}{List of ggplot2 objects (one per \code{tau}).}
  \item{values}{List of data frames with smoothed \code{pred}/\code{obs}
    values per \code{tau}.}
  \item{calib.measures}{List of data frames with ICI, E50, E90, Emax, RSB
    per \code{tau}.}
  \item{OE_summary}{Named numeric vector with \code{OE}, \code{lower},
    \code{upper}.}
}
}
\references{
Van Geloven N et al. (2022). Validation of prediction models in the presence
of competing risks: a guide through modern methods. \emph{BMJ}, 377, e069249.
}
\seealso{
\code{\link{OEComputation}}, \code{\link{score_from_cifs}}
}
