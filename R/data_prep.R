#' @title Data Preparation Utilities
#' @description Functions for preparing and type-coercing data frames for
#'   competing risks benchmarking.
#' @name data_prep
NULL


#' Prepare a data frame for benchmarking
#'
#' Renames the time and event columns to `"time"` and `"event"`, adds a small
#' offset when the minimum observed time is zero (to avoid metric failures),
#' and sorts rows by time.
#'
#' @param df A data frame.
#' @param time_col Name of the column containing observed times.
#' @param event_col Name of the column containing event/status codes
#'   (0 = censored).
#'
#' @return The modified data frame sorted by time.
#' @export
prepare_data <- function(df, time_col, event_col) {
  nm <- names(df)
  nm[nm == time_col]  <- "time"
  nm[nm == event_col] <- "event"
  names(df) <- nm

  min_time <- suppressWarnings(min(df$time, na.rm = TRUE))
  if (is.finite(min_time) && min_time == 0) {
    df$time <- df$time + 0.01
  }

  df[order(df$time), , drop = FALSE]
}


#' Infer column data types from a data frame
#'
#' Returns a schema list compatible with [apply_data_types()].
#'
#' @param df A data frame.
#' @param feature_cols Character vector of feature column names.
#' @param core_cols Character vector of mandatory non-feature columns
#'   (default `c("time","event","row_id")`).
#'
#' @return A named list with elements `version`, `core_cols`,
#'   `feature_cols`, and `types`.
#' @export
infer_data_types <- function(df, feature_cols,
                              core_cols = c("time", "event", "row_id")) {
  types <- stats::setNames(
    vector("list", length(core_cols) + length(feature_cols)),
    c(core_cols, feature_cols)
  )

  types[["time"]]   <- "numeric"
  types[["event"]]  <- "integer"
  types[["row_id"]] <- "character"

  for (cn in feature_cols) {
    if (!cn %in% names(df)) next
    x <- df[[cn]]
    types[[cn]] <- if (is.logical(x))   "logical"
                   else if (is.factor(x))    "factor"
                   else if (is.character(x)) "character"
                   else if (is.integer(x))   "integer"
                   else if (is.numeric(x))   "numeric"
                   else                      "character"
  }

  list(version = 1, core_cols = core_cols,
       feature_cols = feature_cols, types = as.list(types))
}



#' Apply stored data types to a data frame
#'
#' Coerces columns of `df` to the types recorded in `schema` (as produced
#' by [schema_from_metadata()]).
#'
#' @param df A data frame.
#' @param schema A schema list as returned by [schema_from_metadata()].
#'
#' @return The coerced data frame.
#' @export
apply_data_types <- function(df, schema) {
  out   <- df
  types <- schema$types

  if ("time"   %in% names(out)) out$time   <- as.numeric(out$time)
  if ("event"  %in% names(out)) out$event  <- as.integer(out$event)
  if ("row_id" %in% names(out)) out$row_id <- as.character(out$row_id)

  for (cn in names(types)) {
    if (!cn %in% names(out)) next
    if (cn %in% c("time", "event", "row_id")) next
    typ <- types[[cn]]
    out[[cn]] <- switch(typ,
      numeric   = as.numeric(out[[cn]]),
      integer   = as.integer(out[[cn]]),
      logical   = as.logical(out[[cn]]),
      factor    = as.factor(out[[cn]]),
      as.character(out[[cn]])
    )
  }
  out
}
