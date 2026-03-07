\name{prepare_data}
\alias{prepare_data}
\title{Prepare a data frame for benchmarking}
\usage{
prepare_data(df, time_col, event_col)
}
\arguments{
  \item{df}{A data frame.}
  \item{time_col}{Name of the column containing observed times.}
  \item{event_col}{Name of the column containing event/status codes
    (0 = censored).}
}
\description{
Renames time and event columns to \code{"time"} and \code{"event"}, adds a
small positive offset if the minimum observed time is zero (required by some
metric functions), and sorts rows by time.
}
\value{
The modified data frame sorted by \code{time}.
}
