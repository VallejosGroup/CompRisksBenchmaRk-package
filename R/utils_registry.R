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

#' cr_data class
#'
#' An S4 class that holds a validated and pre-processed data frame ready for
#' competing risks model fitting, together with metadata derived during
#' construction.
#'
#' @slot data         A data frame, sorted by \code{time_var} if
#'   \code{sort_by_time = TRUE}.
#' @slot causes       Sorted integer vector of competing cause codes (all
#'   non-zero values found in \code{event_var}).
#' @slot covars       A data frame with columns \code{covars_names} and
#'   \code{covars_types} recording the name and original type of each
#'   feature column.
#' @slot time_var     Name of the time column.
#' @slot event_var    Name of the event/status column.
#' @slot id_var       Name of the subject ID column, or \code{""} if not
#'   supplied.
#' @slot sort_by_time Logical; whether rows were sorted by ascending time.
#' @slot time_offset  Numeric offset that was added to times when the minimum
#'   observed time was zero (0 means no offset was applied).
#'
#' @name cr_data-class
#' @exportClass cr_data
NULL


#' @describeIn cr_data-class Print a concise summary of a \code{cr_data} object.
#' @param object A \code{cr_data} object.
#' @export
methods::setMethod("show", "cr_data", function(object) {
  cat("A cr_data object\n")
  cat(sprintf("  Rows      : %d\n", nrow(object@data)))
  cat(sprintf("  Causes    : %s\n", paste(object@causes, collapse = ", ")))
  cat(sprintf("  time_var  : %s\n", object@time_var))
  cat(sprintf("  event_var : %s\n", object@event_var))
  if (nzchar(object@id_var))
    cat(sprintf("  id_var    : %s\n", object@id_var))
  cat(sprintf("  Covariates: %d  [%s]\n",
              nrow(object@covars),
              paste(object@covars$covars_names, collapse = ", ")))
  if (object@time_offset != 0)
    cat(sprintf("  time_offset applied: %g\n", object@time_offset))
  invisible(object)
})


#' Subset a cr_data object
#'
#' Subsets rows and/or columns of the underlying data frame while preserving
#' all metadata slots.  Causes are re-derived from the row subset so that
#' causes absent in the subset are dropped.  Dropping \code{time_var} or
#' \code{event_var} raises an error; \code{id_var} (if set) is always
#' retained regardless of \code{j}.
#'
#' @param x    A \code{cr_data} object.
#' @param i    Row index (integer, logical, or character vector).  If missing,
#'   all rows are kept.
#' @param j    Column names to keep (character vector).  If missing, all
#'   columns are kept.
#' @param drop Ignored.
#'
#' @return A new \code{cr_data} object.
#' @exportMethod [
methods::setMethod("[", signature("cr_data", "ANY", "ANY", "missing"),
                   function(x, i, j, drop) {
                     new_data <- if (missing(i)) x@data else x@data[i, , drop = FALSE]
                     
                     if (!missing(j)) {
                       protected <- c(x@time_var, x@event_var,
                                      if (nzchar(x@id_var)) x@id_var)
                       bad <- intersect(protected, setdiff(names(new_data), j))
                       if (length(bad))
                         stop(sprintf(
                           "Cannot drop protected column(s): %s",
                           paste(bad, collapse = ", ")
                         ), call. = FALSE)
                       keep_cols <- union(protected, intersect(j, names(new_data)))
                       new_data  <- new_data[, keep_cols, drop = FALSE]
                     }
                     
                     new_covars <- x@covars[x@covars$covars_names %in% names(new_data), ,
                                            drop = FALSE]
                     
                     methods::new("cr_data",
                                  data         = new_data,
                                  causes       = cr_causes(new_data, x@event_var),
                                  covars       = new_covars,
                                  time_var     = x@time_var,
                                  event_var    = x@event_var,
                                  id_var       = x@id_var,
                                  sort_by_time = x@sort_by_time,
                                  time_offset  = x@time_offset
                     )
                   }
)