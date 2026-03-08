#' @title Internal helpers shared across model registrations
#' @name model_helpers
#' @keywords internal
NULL


#' Coerce and validate model inputs
#'
#' Validates `x`, `time`, and `event`, coercing where necessary, and derives
#' the sorted vector of competing causes (all non-zero event codes).
#'
#' @param x Feature matrix or data frame.
#' @param time Numeric vector of observed times.
#' @param event Integer-coercible vector of event codes (0 = censored).
#'
#' @return A named list with elements `x` (data.frame), `time` (numeric),
#'   `event` (integer), and `causes` (sorted integer vector).
#' @noRd
cr_prepare_inputs <- function(x, time, event) {
  if (!is.data.frame(x))
    x <- as.data.frame(x)
  if (!is.numeric(time))
    stop("`time` must be a numeric vector.", call. = FALSE)
  if (!is.integer(event))
    event <- as.integer(as.character(event))
  causes <- sort(unique(event[event != 0]))
  list(x = x, time = time, event = event, causes = causes)
}


#' Build a training data frame and model formula
#'
#' Assembles a data frame from `time`, `event`, and `x`, then builds a
#' model formula with the requested response type.
#'
#' @param x A data frame of features (already coerced).
#' @param time Numeric time vector.
#' @param event Integer event vector.
#' @param response One of `"Hist"` (uses `prodlim::Hist(time, event)`,
#'   default) or `"Surv"` (uses `survival::Surv(time, event)`).
#'
#' @return A named list with elements `dat` (data frame) and `formula`
#'   (the assembled formula object).
#' @noRd
cr_build_formula <- function(x, time, event, response = c("Hist", "Surv")) {
  response <- match.arg(response)
  dat    <- data.frame(time = time, event = event, x, check.names = FALSE)
  covars <- setdiff(names(dat), c("time", "event"))
  f      <- stats::reformulate(covars, response = NULL)
  lhs    <- switch(response,
    Hist = quote(prodlim::Hist(time, event)),
    Surv = quote(survival::Surv(time, event))
  )
  formula <- stats::update(f, paste(deparse(lhs), "~ ."))
  list(dat = dat, formula = formula)
}


#' Extract an argument from an `args` list with a default fallback
#'
#' @param args Named list of arguments (typically passed from the grid).
#' @param name Argument name (character scalar).
#' @param default Value to use when `name` is absent or `NULL` in `args`.
#'
#' @return `args[[name]]` if present and non-`NULL`, otherwise `default`.
#' @noRd
cr_arg <- function(args, name, default = NULL) {
  val <- args[[name]]
  if (is.null(val)) default else val
}
