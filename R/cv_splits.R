#' @title Cross-Validation Split Utilities
#' @description Functions for creating and stratifying nested cross-validation
#'   splits for competing risks benchmarking.
#' @name cv_splits
NULL


#' Create nested cross-validation splits and write them to disk
#'
#' Produces a directory hierarchy suitable for [nested_cv_from_bench()].
#' Outer and inner folds are stratified by event status using
#' [rsample::vfold_cv()].  Fold data are stored as Parquet files; inner
#' fold indices are stored as JSON.
#'
#' @param df A data frame containing at minimum time, event, and covariate
#'   columns.
#' @param time_col Name of the time column.
#' @param event_col Name of the event/status column.
#' @param id_col Optional name of an existing ID column.  If `NULL`, row
#'   numbers are used.
#' @param out_dir Path to the output directory (created if needed).
#' @param outer_folds Number of outer folds (default 5).
#' @param inner_folds Number of inner folds (default 3).
#' @param seed Random seed for reproducibility.
#' @param times Numeric vector of evaluation time points written to
#'   `times.json`.
#'
#' @return Invisibly `TRUE`.
#' @export
create_nested_splits <- function(df,
                                  time_col,
                                  event_col,
                                  id_col       = NULL,
                                  out_dir      = NULL,
                                  outer_folds  = 5,
                                  inner_folds  = 3,
                                  seed         = 123,
                                  times) {
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df <- prepare_data(df = df, time_col = time_col, event_col = event_col)

  if (!is.null(id_col)) {
    df$row_id <- as.character(df[[id_col]])
    df[[id_col]] <- NULL
  } else {
    df$row_id <- as.character(seq_len(nrow(df)))
  }

  core_cols    <- c("time", "event", "row_id")
  feature_cols <- setdiff(names(df), core_cols)

  schema <- infer_data_types(df, feature_cols = feature_cols,
                              core_cols = core_cols)
  jsonlite::write_json(schema,
                       file.path(out_dir, "data_types.json"),
                       auto_unbox = TRUE, pretty = TRUE)

  for_python <- list(
    version     = 1,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    features    = feature_cols,
    time_col    = "time",
    event_col   = "event",
    row_id      = "row_id"
  )
  jsonlite::write_json(for_python,
                       file.path(out_dir, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  jsonlite::write_json(times,
                       file.path(out_dir, "times.json"),
                       auto_unbox = TRUE)

  set.seed(seed)
  outer <- rsample::vfold_cv(df, v = outer_folds, strata = "event")

  for (v in seq_len(outer_folds)) {
    split_v <- outer$splits[[v]]
    train   <- rsample::analysis(split_v)
    test    <- rsample::assessment(split_v)

    out_v_dir <- file.path(out_dir, sprintf("outer_%d", v))
    dir.create(out_v_dir, recursive = TRUE, showWarnings = FALSE)

    arrow::write_parquet(
      dplyr::select(train, dplyr::all_of(c("time", "event", feature_cols, "row_id"))),
      file.path(out_v_dir, "train.parquet")
    )
    arrow::write_parquet(
      dplyr::select(test,  dplyr::all_of(c("time", "event", feature_cols, "row_id"))),
      file.path(out_v_dir, "test.parquet")
    )

    set.seed(seed + v)
    inner     <- rsample::vfold_cv(train, v = inner_folds, strata = "event")
    inner_dir <- file.path(out_v_dir, "inner")
    dir.create(inner_dir, showWarnings = FALSE)

    for (j in seq_len(inner_folds)) {
      split_j    <- inner$splits[[j]]
      train_in   <- rsample::analysis(split_j)
      val_in     <- rsample::assessment(split_j)
      fold_j_dir <- file.path(inner_dir, sprintf("fold_%d", j))
      dir.create(fold_j_dir, showWarnings = FALSE)

      idx_payload <- list(
        train_row_id = train_in[["row_id"]],
        val_row_id   = val_in[["row_id"]]
      )
      jsonlite::write_json(idx_payload,
                           file.path(fold_j_dir, "indices.json"),
                           auto_unbox = TRUE)
    }
  }

  invisible(TRUE)
}


#' Stratified k-fold CV split object
#'
#' Creates an `rsample` v-fold cross-validation object stratified by the
#' event column so that each fold contains a representative mix of event
#' types.
#'
#' @param df A data frame containing an `event` column.
#' @param time  (Unused; kept for signature compatibility.)
#' @param event (Unused; `df$event` is used directly.)
#' @param folds Number of folds (default 5).
#'
#' @return An `rsample::vfold_cv` object.
#' @export
stratification <- function(df, time, event, folds = 5) {
  df$event <- as.factor(df$event)
  rsample::vfold_cv(df, v = folds, strata = "event")
}
