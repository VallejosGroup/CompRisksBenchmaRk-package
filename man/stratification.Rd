\name{stratification}
\alias{stratification}
\title{Stratified k-fold cross-validation splits}
\usage{
stratification(df, time, event, folds = 5)
}
\arguments{
  \item{df}{A data frame containing an \code{event} column.}
  \item{time}{Unused; retained for signature compatibility.}
  \item{event}{Unused; \code{df$event} is used directly.}
  \item{folds}{Number of folds (default 5).}
}
\description{
Creates an \code{rsample} v-fold cross-validation object stratified by the
event column so that each fold contains a representative proportion of each
event type.
}
\value{
An \code{rsample::vfold_cv} object.
}
