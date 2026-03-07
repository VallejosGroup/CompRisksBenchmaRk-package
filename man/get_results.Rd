\name{get_results}
\alias{get_results}
\title{Load and score benchmark predictions across models}
\usage{
get_results(
  data_root,
  dataset,
  dataset_prior_imputation = NULL,
  python_models = NULL,
  r_models,
  outer_folds,
  horizon,
  seed,
  cause_of_interest,
  sample_curves
)
}
\arguments{
  \item{data_root}{Root directory containing dataset subdirectories.}
  \item{dataset}{Name of the dataset subdirectory.}
  \item{dataset_prior_imputation}{Name of the directory holding the original
    (pre-imputation) splits, used to locate the test sets.  Required when
    \code{dataset} contains the string \code{"imputed"}.}
  \item{python_models}{Character vector of model keys whose CIF arrays are
    stored in \code{[K, Tm, n]} order (Python convention).}
  \item{r_models}{Character vector of model keys whose CIF arrays are in
    \code{[n, K, Tm]} order.}
  \item{outer_folds}{Number of outer CV folds.}
  \item{horizon}{Scalar evaluation time used for calibration and O/E plots.}
  \item{seed}{Random seed for reproducible curve sampling.}
  \item{cause_of_interest}{Integer index of the cause to evaluate (1-based).}
  \item{sample_curves}{Logical; if \code{TRUE} a sample of CIF trajectories
    is printed for each fold.}
}
\description{
Reads pre-computed CIF Parquet files from
\file{data_root/dataset/<model>/outer_v/cif.parquet}, scores them with
\code{\link{score_from_cifs}}, plots calibration curves via
\code{\link{CalibrationPlot}}, and returns all results and predictions
in a single list.
}
\value{
A named list with elements:
\describe{
  \item{all_test_folds}{Data frame of all outer test observations concatenated
    across folds, with extra columns for each model's predicted CIF at
    \code{horizon}.}
  \item{predictor_cols}{Character vector naming the CIF columns added to
    \code{all_test_folds}.}
  \item{all_cifs}{Named list (by model) of CIF arrays concatenated across
    folds.}
  \item{all_rmlt}{Named list (by model) of restricted mean lifetime matrices.}
  \item{res}{Named list (by model) of per-fold scoring results.}
  \item{res_rmlt}{Named list (by model) of per-fold RMLT C-index results.}
  \item{oes}{Named list of O/E summaries per model at \code{horizon}.}
}
}
\seealso{
\code{\link{nested_cv_from_bench}}, \code{\link{score_from_cifs}},
\code{\link{CalibrationPlot}}
}
