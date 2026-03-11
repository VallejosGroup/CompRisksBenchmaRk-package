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
#' @slot cens_code    Integer code used to denote censored observations in
#'   \code{event_var} (default \code{0L}).
#'
#' @name cr_data-class
#' @exportClass cr_data
NULL

#' @noRd
methods::setClass("cr_data", representation(
  data         = "data.frame",
  causes       = "integer",
  covars       = "data.frame",
  time_var     = "character",
  event_var    = "character",
  id_var       = "character",
  sort_by_time = "logical",
  time_offset  = "numeric",
  cens_code    = "integer"
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
#' @param cens_code  Integer code denoting censored observations in
#'   \code{event_var} (default \code{0}).
#'
#' @return A \code{cr_data} S4 object.
#' @export
cr_data <- function(data, time_var, event_var,
                    sort_by_time = TRUE,
                    id_var       = NULL,
                    time_offset  = 0,
                    cens_code    = 0L) {

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

  if (!is.numeric(cens_code) || length(cens_code) != 1)
    stop("`cens_code` must be a single integer value.", call. = FALSE)
  cens_code <- as.integer(cens_code)

  # --- Zero-time guard ---
  min_time <- suppressWarnings(min(data[[time_var]], na.rm = TRUE))
  if (is.finite(min_time) && min_time == 0) {
    if (time_offset == 0)
      stop(
        "Minimum observed time is zero, which can cause metric failures. ",
        "Either exclude observations with ", time_var, " = 0, or re-run ",
        "with a positive `time_offset` (e.g. time_offset = 0.01).",
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
  causes <- .cr_causes(data, event_var, cens_code)
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
    time_offset  = time_offset,
    cens_code    = cens_code
  )
}


#' Summarise a cr_data object
#'
#' Produces an HTML summary table of all covariates and the observed survival
#' times, stratified by event status (censored + each competing cause), using
#' the \pkg{table1} package.
#'
#' @param object A \code{cr_data} object.
#' @param ... Currently unused.
#'
#' @return A \code{table1} object (printed as HTML, invisibly returned).
#' @export
methods::setMethod("summary", "cr_data", function(object, ...) {
  d         <- object@data
  time_var  <- object@time_var
  event_var <- object@event_var
  causes    <- object@causes
  cens_code <- object@cens_code

  # Build labelled event factor: "Censored" + "Cause k" for each cause
  ev_int    <- as.integer(d[[event_var]])
  cens_lbl  <- paste0("Censored (", cens_code, ")")
  cause_lbls <- paste0("Cause ", causes)
  lvls      <- c(cens_code, causes)
  lbls      <- c(cens_lbl, cause_lbls)
  ev_factor <- factor(ev_int, levels = lvls, labels = lbls)

  # Assemble table data: time + covariates + event factor
  covar_names <- object@covars$covars_names
  tbl_data    <- cbind(
    d[, c(time_var, covar_names), drop = FALSE],
    .event_status = ev_factor
  )

  # Apply table1 labels: time_var gets a readable label
  table1::label(tbl_data[[time_var]]) <- paste0("Survival time (", time_var, ")")

  # Build formula: all vars on LHS, event factor on RHS
  rhs_vars <- c(time_var, covar_names)
  f <- stats::as.formula(
    paste("~", paste(rhs_vars, collapse = " + "), "| .event_status")
  )

  tbl <- table1::table1(f, data = tbl_data, overall = "Overall")
  print(tbl)
  invisible(tbl)
})


#' @describeIn cr_data-class Print a concise summary of a \code{cr_data} object.
#' @param object A \code{cr_data} object.
#' @export
methods::setMethod("show", "cr_data", function(object) {
  cat("A cr_data object\n")
  cat(sprintf("  Rows      : %d\n", nrow(object@data)))
  cat(sprintf("  Columns   : %d\n", ncol(object@data)))
  cat(sprintf("  Causes    : %s\n", paste(object@causes, collapse = ", ")))
  cat(sprintf("  cens_code : %d\n", object@cens_code))
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


#' @describeIn cr_data-class Dimensions of the underlying data frame.
#' @param x A \code{cr_data} object.
#' @export
methods::setMethod("dim", "cr_data", function(x) dim(x@data))


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
      causes       = .cr_causes(new_data, x@event_var, x@cens_code),
      covars       = new_covars,
      time_var     = x@time_var,
      event_var    = x@event_var,
      id_var       = x@id_var,
      sort_by_time = x@sort_by_time,
      time_offset  = x@time_offset,
      cens_code    = x@cens_code
    )
  }
)
