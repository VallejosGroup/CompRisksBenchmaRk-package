\name{impute_split}
\alias{impute_split}
\title{Impute missing values in a train/test split}
\usage{
impute_split(
  train,
  test,
  feature_cols,
  method = c("none", "median_mode", "mice"),
  m = 1L,
  maxit = 100L,
  seed = 1L,
  drop_cols_after_impute = NULL
)
}
\arguments{
  \item{train}{Training data frame.}
  \item{test}{Test (or validation) data frame.}
  \item{feature_cols}{Character vector of feature column names to impute.}
  \item{method}{One of \code{"none"} (passthrough), \code{"median_mode"}
    (single imputation), or \code{"mice"} (multiple imputation).}
  \item{m}{Number of multiply-imputed datasets (only used for
    \code{method = "mice"}).}
  \item{maxit}{Number of MICE iterations.}
  \item{seed}{Integer seed for MICE.}
  \item{drop_cols_after_impute}{Optional character vector of column names to
    drop after imputation.}
}
\description{
Detects missing values in \code{feature_cols} and applies the chosen
imputation strategy.  Imputation parameters are learnt from \code{train} and
applied to both \code{train} and \code{test} to avoid data leakage.
}
\value{
A list with:
\describe{
  \item{m}{Number of imputed datasets.}
  \item{imputations}{A list of length \code{m}, each a list with
    \code{$train} and \code{$test} data frames.}
}
}
\seealso{
\code{\link{impute_mice}}
}
