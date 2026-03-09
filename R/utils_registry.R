#' @title Internal registry utilities
#' @description Internal environment storing registered models, and helper
#'   functions shared across model registrations.
#' @name utils_registry
#' @keywords internal
NULL

# Internal environment that stores all registered models.
.cr_models <- new.env(parent = emptyenv())


#' Extract competing causes from a data frame
#'
#' Derives the sorted vector of competing causes (all non-zero event codes)
#' from the named event column of a data frame, coercing to integer if needed.
#'
#' @param data A data frame containing the event column.
#' @param event_var Character string naming the event column in `data`.
#'
#' @return A sorted integer vector of cause codes.
#' @noRd
cr_causes <- function(data, event_var) {
  event <- data[[event_var]]
  if (!is.integer(event))
    event <- as.integer(as.character(event))
  sort(unique(event[event != 0]))
}


#' @noRd
methods::setClass("cr_data", representation(
  data         = "data.frame",
  causes       = "integer",
  covars       = "data.frame",
  time_var     = "character",
  event_var    = "character",
  id_var       = "character",
  sort_by_time = "logical",
  time_offset  = "numeric"
))


#' Prepare data for competing risks model fitting
#'
#' Validates inputs, coerces the event column to integer, records covariate
#' types, and optionally sorts by time.  Returns a \code{cr_data} object
#' consumed by model \code{fit} functions.
#'
#' @param data        A data frame.
#' @param time_var    Name of the time column (must be numeric, no NAs).
#' @param event_var   Name of the event/status column (numeric, integer, or
#'   factor; 0 = censored; no NAs).
#' @param sort_by_time Logical; sort rows by ascending time (default \code{TRUE}).
#' @param id_var      Optional name of the subject ID column.  Not coerced.
#' @param time_offset Non-negative numeric offset added to all times when the
#'   minimum observed time is zero.  Must be \code{> 0} in that case; an
#'   error is raised when it is \code{0} (the default).
#'
#' @return A \code{cr_data} S4 object.
#' @export
cr_data <- function(data, time_var, event_var,
                    sort_by_time = TRUE,
                    id_var       = NULL,
                    time_offset  = 0) {

  # --- Input validation ---
  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  if (nrow(data) == 0)
    stop("`data` has no rows.", call. = FALSE)

  if (!is.character(time_var) || length(time_var) != 1 || !nzchar(time_var))
    stop("`time_var` must be a single non-empty string.", call. = FALSE)
  if (!is.character(event_var) || length(event_var) != 1 || !nzchar(event_var))
    stop("`event_var` must be a single non-empty string.", call. = FALSE)

  if (!time_var %in% names(data))
    stop(sprintf("`time_var` '%s' not found in `data`.", time_var), call. = FALSE)
  if (!event_var %in% names(data))
    stop(sprintf("`event_var` '%s' not found in `data`.", event_var), call. = FALSE)

  if (!is.null(id_var)) {
    if (!is.character(id_var) || length(id_var) != 1 || !nzchar(id_var))
      stop("`id_var` must be a single non-empty string.", call. = FALSE)
    if (!id_var %in% names(data))
      stop(sprintf("`id_var` '%s' not found in `data`.", id_var), call. = FALSE)
    if (any(is.na(data[[id_var]])))
      stop(sprintf("`id_var` '%s' contains NA values.", id_var), call. = FALSE)
  }

  if (!is.numeric(data[[time_var]]))
    stop(sprintf("`time_var` '%s' must be numeric.", time_var), call. = FALSE)
  if (any(is.na(data[[time_var]])))
    stop(sprintf("`time_var` '%s' contains NA values.", time_var), call. = FALSE)

  event <- data[[event_var]]
  if (!is.numeric(event) && !is.integer(event) && !is.factor(event))
    stop(sprintf("`event_var` '%s' must be numeric, integer, or factor.", event_var),
         call. = FALSE)
  if (any(is.na(event)))
    stop(sprintf("`event_var` '%s' contains NA values.", event_var), call. = FALSE)

  if (!is.numeric(time_offset) || length(time_offset) != 1 || time_offset < 0)
    stop("`time_offset` must be a single non-negative number.", call. = FALSE)

  # --- Zero-time guard ---
  min_time <- suppressWarnings(min(data[[time_var]], na.rm = TRUE))
  if (is.finite(min_time) && min_time == 0) {
    if (time_offset == 0)
      stop(
        "Minimum observed time is zero, which can cause metric failures. ",
        "Re-run with a positive `time_offset` (e.g. time_offset = 0.01).",
        call. = FALSE
      )
    data[[time_var]] <- data[[time_var]] + time_offset
  }

  # --- Type coercion ---
  data[[event_var]] <- as.integer(as.character(data[[event_var]]))

  feature_cols <- setdiff(names(data), c(time_var, event_var, id_var))
  covars_types <- vector("character", length(feature_cols))
  for (i in seq_along(feature_cols)) {
    x <- data[[feature_cols[[i]]]]
    if (is.logical(x)) {
      covars_types[[i]] <- "logical"
    } else if (is.factor(x)) {
      covars_types[[i]] <- "factor"
    } else if (is.character(x)) {
      covars_types[[i]] <- "character"
    } else if (is.integer(x)) {
      covars_types[[i]] <- "integer"
    } else {
      covars_types[[i]]         <- "numeric"
      data[[feature_cols[[i]]]] <- as.numeric(x)
    }
  }

  # --- Optional sort by time ---
  if (sort_by_time)
    data <- data[order(data[[time_var]]), , drop = FALSE]

  # --- Cause and covariate extraction ---
  causes <- cr_causes(data, event_var)
  covars <- data.frame(covars_names = feature_cols,
                       covars_types = covars_types,
                       stringsAsFactors = FALSE)

  methods::new("cr_data",
    data         = data,
    causes       = causes,
    covars       = covars,
    time_var     = time_var,
    event_var    = event_var,
    id_var       = if (!is.null(id_var)) id_var else "",
    sort_by_time = sort_by_time,
    time_offset  = time_offset
  )
}
