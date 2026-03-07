\name{register_cr_model}
\alias{register_cr_model}
\title{Register a competing risks model}
\usage{
register_cr_model(key, fit, predict_cif, info)
}
\arguments{
  \item{key}{Character string identifying the model (e.g. \code{"FGR"}).}
  \item{fit}{A function \code{function(x, y_time, y_event, args)} that trains
    the model and returns a fitted object.}
  \item{predict_cif}{A function \code{function(fit_obj, x_new, times)} that
    returns an array of dimensions \code{[n, K, Tm]} of predicted CIFs.}
  \item{info}{A zero-argument function returning a named list with at least
    \code{name}, \code{supports}, \code{needs_tuning}, and
    \code{default_grid}.}
}
\description{
Stores a competing risks model specification (fit, predict, info) in the
package-internal registry environment so it can be retrieved by
\code{\link{get_cr_model}} and used inside \code{\link{nested_cv_from_bench}}.
If a model with the same \code{key} already exists a warning is issued and
the entry is overwritten.
}
\value{
Invisibly returns \code{NULL}. Called for side effects.
}
\seealso{
\code{\link{get_cr_model}}, \code{\link{list_cr_models}}
}
