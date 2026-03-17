#' Create nested cross-validation splits and optionally write them to disk
#'
#' Produces outer and inner folds stratified by event status using
#' [rsample::vfold_cv()].  Always returns a list of splits; optionally also
#' writes fold data as Parquet files and inner indices as JSON to a directory
#' hierarchy suitable for [nested_cv_from_bench()].
#'
#' @param cr          A [cr_data()] object.
#' @param store_dir   Path to the output directory (created if needed).
#'   Required when `store = TRUE`.
#' @param outer_folds Number of outer folds (default `5`).
#' @param inner_folds Number of inner folds (default `3`).
#' @param seed        Random seed for reproducibility.
#' @param store       Logical; if `TRUE` (default `FALSE`), write splits and
#'   manifest to `store_dir` in addition to returning the list.
#'
#' @return A named list with one element per outer fold (named `outer_1`,
#'   `outer_2`, …). Each element is a list with:
#'   \describe{
#'     \item{`train`}{Training data frame.}
#'     \item{`test`}{Test data frame.}
#'     \item{`inner`}{A named list with one element per inner fold
#'       (`inner_1`, …), each containing `train_ids` and `val_ids` — character
#'       vectors of `id_var` values for the inner training and validation sets.}
#'   }
#' @export
create_nested_splits <- function(cr,
                                 store_dir    = NULL,
                                 outer_folds  = 5,
                                 inner_folds  = 3,
                                 seed         = 123,
                                 store        = FALSE) {
  .check_cr(cr)
  if (store) {
    if (is.null(store_dir) || !nzchar(store_dir))
      stop("`store_dir` must be a non-empty character string when `store = TRUE`.",
           call. = FALSE)
    if (!dir.exists(store_dir))
      dir.create(store_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  time_var     <- cr@time_var
  event_var    <- cr@event_var
  id_var       <- cr@id_var
  feature_vars <- cr@feature_vars
  
  df <- cr@data
  
  set.seed(seed)
  outer_cv <- rsample::vfold_cv(df, v = outer_folds, strata = event_var)
  
  result <- stats::setNames(vector("list", outer_folds),
                            sprintf("outer_%d", seq_len(outer_folds)))
  
  for (v in seq_len(outer_folds)) {
    split_v <- outer_cv$splits[[v]]
    train   <- rsample::analysis(split_v)
    test    <- rsample::assessment(split_v)
    
    set.seed(seed + v)
    inner_cv <- rsample::vfold_cv(train, v = inner_folds, strata = event_var)
    
    inner_list <- stats::setNames(vector("list", inner_folds),
                                  sprintf("inner_%d", seq_len(inner_folds)))
    for (j in seq_len(inner_folds)) {
      split_j  <- inner_cv$splits[[j]]
      train_in <- rsample::analysis(split_j)
      val_in   <- rsample::assessment(split_j)
      inner_list[[j]] <- list(
        train_ids = train_in[[id_var]],
        val_ids   = val_in[[id_var]]
      )
    }
    
    result[[v]] <- list(train = train, test = test, inner = inner_list)
    
    if (store) {
      out_v_dir <- file.path(store_dir, sprintf("outer_%d", v))
      dir.create(out_v_dir, recursive = TRUE, showWarnings = FALSE)
      
      arrow::write_parquet(train, file.path(out_v_dir, "train.parquet"))
      arrow::write_parquet(test,  file.path(out_v_dir, "test.parquet"))
      
      inner_dir <- file.path(out_v_dir, "inner")
      dir.create(inner_dir, showWarnings = FALSE)
      
      for (j in seq_len(inner_folds)) {
        fold_j_dir <- file.path(inner_dir, sprintf("fold_%d", j))
        dir.create(fold_j_dir, showWarnings = FALSE)
        jsonlite::write_json(
          list(train_row_id = inner_list[[j]]$train_ids,
               val_row_id   = inner_list[[j]]$val_ids),
          file.path(fold_j_dir, "indices.json"),
          auto_unbox = TRUE
        )
      }
    }
  }
  
  if (store) {
    for_python <- list(
      version     = 1,
      outer_folds = outer_folds,
      inner_folds = inner_folds,
      features    = feature_vars,
      time_col    = time_var,
      event_col   = event_var,
      row_id      = id_var
    )
    jsonlite::write_json(for_python,
                         file.path(store_dir, "manifest.json"),
                         auto_unbox = TRUE, pretty = TRUE)
  }
  
  result
}


#' Convert cr_metadata to legacy schema format
#'
#' Converts the metadata list returned by [cr_metadata()] into the schema
#' format previously produced by `infer_data_types()`, for backwards
#' compatibility with pipelines that read `data_types.json`.
#'
#' @param metadata A named list as returned by [cr_metadata()].
#'
#' @return A named list with elements `version`, `core_cols`,
#'   `feature_cols`, and `types`, compatible with the legacy
#'   `data_types.json` format.
#' @export
schema_from_metadata <- function(metadata) {
  core_cols    <- c(metadata$time_var, metadata$event_var, metadata$id_var)
  feature_vars <- metadata$feature_vars
  feature_types <- metadata$feature_types
  
  types <- as.list(stats::setNames(
    c(rep(NA_character_, length(core_cols)), feature_types),
    c(core_cols, feature_vars)
  ))
  
  list(
    version      = 1,
    core_cols    = core_cols,
    feature_cols = feature_vars,
    types        = types
  )
}