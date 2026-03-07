\name{nested_cv_from_bench}
\alias{nested_cv_from_bench}
\title{Nested cross-validation using pre-defined benchmark splits}
\usage{
nested_cv_from_bench(
  out_dir = "../BenchResults",
  model_key,
  grid = NULL,
  seed = 123,
  verbose = TRUE,
  drop_cols_after_impute = NULL
)
}
\arguments{
  \item{out_dir}{Path to the benchmark root directory containing
    \file{manifest.json}, \file{times.json}, \file{data_types.json}, and
    \file{outer_*/} fold subdirectories.}
  \item{model_key}{Character string key of a registered model
    (see \code{\link{register_cr_model}}).}
  \item{grid}{Named list of hyperparameter values to search, or \code{NULL}
    for models that need no tuning. Expanded via \code{tidyr::expand_grid}.}
  \item{seed}{Integer random seed (added to the fold index for
    reproducibility).}
  \item{verbose}{Logical; print progress messages.}
  \item{drop_cols_after_impute}{Optional character vector of feature column
    names to remove after imputation is complete.}
}
\description{
Runs a full nested cross-validation workflow against a pre-built benchmark
directory layout produced by the companion Python pipeline.  For each outer
fold the function:
\enumerate{
  \item Loads (possibly multiply-imputed) train/test parquet files.
  \item If the model has \code{needs_tuning = TRUE}, performs inner-fold
    hyperparameter selection minimising the integrated Brier score per cause.
  \item Refits the model on the full outer training set with the best
    hyperparameters and predicts CIFs on the outer test set.
  \item Saves the CIF predictions as \file{<out_dir>/<model_key>/outer_v/cif.parquet}.
}
}
\value{
A list of length \code{outer_folds}.  Each element contains \code{fold},
\code{model_key}, \code{times}, and \code{metrics} (including per-cause IBS).
For models with tuning, \code{best_cfg} is also included.
}
\seealso{
\code{\link{get_results}}, \code{\link{score_from_cifs}}
}
