#' Create nested cross-validation splits and write them to disk
#'
#' Produces a directory hierarchy suitable for [nested_cv_from_bench()].
#' Outer and inner folds are stratified by event status using
#' [rsample::vfold_cv()].  Fold data are stored as Parquet files; inner
#' fold indices are stored as JSON.
#'
#' @param cr          A [cr_data()] object.
#' @param out_dir     Path to the output directory (created if needed).
#' @param outer_folds Number of outer folds (default `5`).
#' @param inner_folds Number of inner folds (default `3`).
#' @param seed        Random seed for reproducibility.
#' @param times       Numeric vector of evaluation time points written to
#'   `times.json`.
#'
#' @return Invisibly `TRUE`.
#' @export
create_nested_splits <- function(cr,
                                  out_dir      = NULL,
                                  outer_folds  = 5,
                                  inner_folds  = 3,
                                  seed         = 123,
                                  times) {
  .check_cr(cr)
  if (is.null(out_dir) || !nzchar(out_dir))
    stop("`out_dir` must be a non-empty character string.", call. = FALSE)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  time_var    <- cr@time_var
  event_var   <- cr@event_var
  id_var      <- cr@id_var
  feature_vars <- cr@feature_vars

  # Build a flat data frame with standardised column names for the pipeline
  df          <- as.data.frame(cr@data)
  df[["time"]]   <- df[[time_var]]
  df[["event"]]  <- df[[event_var]]
  df[["row_id"]] <- if (!is.null(id_var) && id_var %in% names(df))
                      as.character(df[[id_var]])
                    else
                      as.character(seq_len(nrow(df)))

  # Add small offset if min time is zero to avoid metric failures
  min_time <- suppressWarnings(min(df[["time"]], na.rm = TRUE))
  if (is.finite(min_time) && min_time == 0)
    df[["time"]] <- df[["time"]] + 0.01

  df <- df[order(df[["time"]]), , drop = FALSE]

  core_cols <- c("time", "event", "row_id")

  schema <- infer_data_types(df, feature_cols = feature_vars,
                              core_cols = core_cols)
  jsonlite::write_json(schema,
                       file.path(out_dir, "data_types.json"),
                       auto_unbox = TRUE, pretty = TRUE)

  for_python <- list(
    version     = 1,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    features    = feature_vars,
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

    keep_cols <- c("time", "event", feature_vars, "row_id")
    arrow::write_parquet(dplyr::select(train, dplyr::all_of(keep_cols)),
                         file.path(out_v_dir, "train.parquet"))
    arrow::write_parquet(dplyr::select(test,  dplyr::all_of(keep_cols)),
                         file.path(out_v_dir, "test.parquet"))

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

      jsonlite::write_json(
        list(train_row_id = train_in[["row_id"]],
             val_row_id   = val_in[["row_id"]]),
        file.path(fold_j_dir, "indices.json"),
        auto_unbox = TRUE
      )
    }
  }

  invisible(TRUE)
}
