#' @title Internal helpers shared across model registrations
#' @name model_helpers
#' @keywords internal
NULL


#' Coerce and validate model inputs
#'
#' Converts `x`, `y_time`, and `y_event` to the expected types and derives
#' the sorted vector of competing causes (all non-zero event codes).
#'
#' @param x Feature matrix or data frame.
#' @param y_time Numeric vector of observed times.
#' @param y_event Integer-coercible vector of event codes (0 = censored).
#'
#' @return A named list with elements `x` (data.frame), `y_time` (numeric),
#'   `y_event` (integer), and `causes` (sorted integer vector).
#' @noRd
cr_prepare_inputs <- function(x, y_time, y_event) {
  x       <- as.data.frame(x)
  y_time  <- as.numeric(y_time)
  y_event <- as.integer(as.character(y_event))
  causes  <- sort(unique(y_event[y_event != 0L]))
  list(x = x, y_time = y_time, y_event = y_event, causes = causes)
}


#' Build a training data frame and model formula
#'
#' Assembles a data frame from `y_time`, `y_event`, and `x`, then builds a
#' model formula with the requested response type.
#'
#' @param x A data frame of features (already coerced).
#' @param y_time Numeric time vector.
#' @param y_event Integer event vector.
#' @param response One of `"Hist"` (uses `prodlim::Hist(time, event)`,
#'   default) or `"Surv"` (uses `survival::Surv(time, event)`).
#'
#' @return A named list with elements `dat` (data frame) and `formula`
#'   (the assembled formula object).
#' @noRd
cr_build_formula <- function(x, y_time, y_event, response = c("Hist", "Surv")) {
  response <- match.arg(response)
  dat    <- data.frame(time = y_time, event = y_event, x, check.names = FALSE)
  covars <- setdiff(names(dat), c("time", "event"))
  f      <- stats::reformulate(covars, response = NULL)
  lhs    <- switch(response,
                   Hist = quote(prodlim::Hist(time, event)),
                   Surv = quote(survival::Surv(time, event))
  )
  formula <- stats::update(f, paste(deparse(lhs), "~ ."))
  list(dat = dat, formula = formula)
}


#' Initialise the output CIF array and coerce prediction inputs
#'
#' @param fit_obj A fitted model object with an element `causes`.
#' @param x_new Feature matrix or data frame for new observations.
#' @param times Numeric vector of evaluation times.
#'
#' @return A named list with elements `out` (zero-filled array of dim
#'   `[n, K, Tm]`), `x_new` (data.frame), `times` (numeric), `n`, and `K`.
#' @noRd
cr_init_cif_array <- function(fit_obj, x_new, times) {
  x_new <- as.data.frame(x_new)
  times <- as.numeric(times)
  n     <- nrow(x_new)
  K     <- length(fit_obj$causes)
  out   <- array(0, dim = c(n, K, length(times)))
  list(out = out, x_new = x_new, times = times, n = n, K = K)
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