\name{sim_cmprks}
\alias{sim_cmprks}
\title{Simulate a competing risks dataset}
\usage{
sim_cmprks(n = 500, seed = 1)
}
\arguments{
  \item{n}{Number of subjects to simulate.}
  \item{seed}{Random seed.}
}
\description{
Generates a simple synthetic competing risks dataset with two causes and
random covariates, useful for testing and vignettes.
}
\value{
A data frame with columns \code{time}, \code{event}, and covariate columns.
}
