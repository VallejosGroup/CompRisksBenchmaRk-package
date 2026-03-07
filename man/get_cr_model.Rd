\name{get_cr_model}
\alias{get_cr_model}
\title{Retrieve a registered competing risks model}
\usage{
get_cr_model(key)
}
\arguments{
  \item{key}{Character string model key.}
}
\description{
Looks up a model previously registered with \code{\link{register_cr_model}}
and returns its specification list.
}
\value{
A list with elements \code{fit}, \code{predict_cif}, and \code{info}.
}
\seealso{
\code{\link{register_cr_model}}, \code{\link{list_cr_models}}
}
