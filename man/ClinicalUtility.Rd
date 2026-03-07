\name{ClinicalUtility}
\alias{ClinicalUtility}
\title{Decision curve analysis for competing risks}
\usage{
ClinicalUtility(
  data,
  outcome,
  ttoutcome,
  timepoint,
  predictors,
  xstart = 0.01,
  xstop = 0.99,
  xby = 0.01,
  ymin = -0.05,
  probability = NULL,
  harm = NULL,
  graph = TRUE,
  intervention = FALSE,
  interventionper = 100,
  smooth = FALSE,
  loess.span = 0.10,
  cmprsk = TRUE
)
}
\arguments{
  \item{data}{Data frame of complete cases.}
  \item{outcome}{Name of the event/status column.}
  \item{ttoutcome}{Name of the time-to-outcome column.}
  \item{timepoint}{Scalar evaluation time for the decision curve.}
  \item{predictors}{Character vector of predicted risk column names
    (must lie in \code{[0, 1]}).}
  \item{xstart}{Lower threshold probability (default 0.01).}
  \item{xstop}{Upper threshold probability (default 0.99).}
  \item{xby}{Step size for threshold grid (default 0.01).}
  \item{ymin}{Minimum y-axis value for plots.}
  \item{probability}{Logical vector indicating which predictors are already
    on the probability scale.}
  \item{harm}{Numeric vector of model-specific harm values.}
  \item{graph}{Logical; produce a ggplot2 plot.}
  \item{intervention}{Logical; plot interventions avoided instead of net
    benefit.}
  \item{interventionper}{Denominator for interventions avoided
    (default 100).}
  \item{smooth}{Logical; apply LOESS smoothing to the curves.}
  \item{loess.span}{LOESS span parameter (default 0.10).}
  \item{cmprsk}{Logical; use competing-risks cumulative incidence
    (\code{cmprsk::cuminc}) instead of Kaplan-Meier.}
}
\description{
Implements decision curve analysis (DCA) adapted for competing risks,
following Giardiello et al. / Van Geloven et al. (2022).  Returns net
benefit and interventions-avoided curves together with a ggplot2 figure.
}
\value{
A list with \code{N}, \code{predictors}, \code{net.benefit},
\code{interventions.avoided}, and \code{plot} (a ggplot2 object when
\code{graph = TRUE}).
}
\references{
Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method for
evaluating prediction models. \emph{Med Decis Making}, 26(6):565-574.
}
\seealso{
\code{\link{plot_dca_gg}}
}
