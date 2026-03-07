\name{CompRisksBenchmaRk-package}
\alias{CompRisksBenchmaRk-package}
\alias{CompRisksBenchmaRk}
\docType{package}
\title{Competing Risks Benchmarking Framework}
\description{
A benchmarking framework for competing risks survival models supporting
nested cross-validation, multiple imputation, and a suite of evaluation
metrics (Brier score, AUC, C-index, calibration, DCA).

Built-in models: Fine-Gray (\code{"FGR"}), penalised Fine-Gray
(\code{"FGRP"}), cause-specific Cox PH (\code{"csCPH"}), random survival
forest (\code{"RSF"}), and Vertical Model (\code{"VM_spline"},
\code{"VM_piecewise"}).  Additional models can be added via
\code{\link{register_cr_model}}.
}
\author{
VallejosGroup
}
\references{
Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution
of a competing risk. \emph{JASA}, 94(446):496-509.

Nicolaie MA et al. (2010). Dynamic prediction by landmarking in competing
risks. \emph{Statistics in Medicine}, 29(12):1190-1205.

Van Geloven N et al. (2022). Validation of prediction models in the presence
of competing risks. \emph{BMJ}, 377, e069249.
}
\keyword{package}
