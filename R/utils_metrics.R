#' @title Internal metric utilities
#' @description Internal helper functions shared across metric computations.
#' @name utils_metrics
#' @keywords internal
NULL


#' Check that `cr` is a cr_data object
#'
#' @param cr Object to check.
#' @noRd
.check_cr <- function(cr) {
  if (!methods::is(cr, "cr_data"))
    stop("`cr` must be a cr_data object.", call. = FALSE)
  invisible(NULL)
}


#' Resolve and validate the `cif` / `fit` / `cif_time_grid` arguments
#'
#' Enforces mutual exclusion of `cif` and `fit`, checks that `cif_time_grid`
#' is `NULL` when `cif` is supplied directly, and calls [predict_cif()] when
#' `fit` is supplied.  Returns the resolved `cif` list.
#'
#' @param cif           A predict_cif() list or `NULL`.
#' @param fit           A fitted model object or `NULL`.
#' @param cif_time_grid Numeric time grid (required when `fit` is non-`NULL`).
#' @param cr            A `cr_data` object (passed to [predict_cif()]).
#'
#' @return The resolved `cif` list (as returned by [predict_cif()]).
#' @noRd
.resolve_cif <- function(cif, fit, cif_time_grid, cr) {
  if (is.null(cif) && is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are non-NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(cif_time_grid))
    stop("`cif_time_grid` must be NULL when `cif` is supplied directly.",
         call. = FALSE)
  if (!is.null(fit)) {
    if (is.null(cif_time_grid))
      stop("`cif_time_grid` must be provided when `fit` is non-NULL.",
           call. = FALSE)
    cif <- predict_cif(fit, newdata = cr, time_grid = cif_time_grid)
  }
  cif
}


#' Validate a cif list and unpack into array + time grid
#'
#' Checks that `cif` is a list with a 3-D numeric `$cif` array and a
#' non-empty numeric `$time_grid`, and that the number of cause dimensions
#' matches `cr@causes`.  Returns a list with elements `cif` (the array) and
#' `time_grid` (the numeric vector).
#'
#' @param cif A list as returned by [predict_cif()].
#' @param cr  A `cr_data` object.
#'
#' @return A list with elements `cif` (3-D array), `time_grid` (numeric
#'   vector), and `model_key` (character, `NULL` if not present in input).
#' @noRd
.validate_and_unpack_cif <- function(cif, cr) {
  if (!is.list(cif) || !all(c("cif", "time_grid") %in% names(cif)))
    stop("`cif` must be a list with elements `$cif` and `$time_grid`, ",
         "as returned by predict_cif().", call. = FALSE)
  if (!is.numeric(cif$cif) || length(dim(cif$cif)) != 3)
    stop("`cif$cif` must be a 3-D numeric array `[n, K, Tm]`.", call. = FALSE)
  if (!is.numeric(cif$time_grid) || length(cif$time_grid) == 0)
    stop("`cif$time_grid` must be a non-empty numeric vector.", call. = FALSE)
  if (dim(cif$cif)[2] != length(cr@causes))
    stop(sprintf(
      "`cif` has %d cause dimension(s) but `cr` contains %d cause(s) (%s).",
      dim(cif$cif)[2], length(cr@causes), paste(cr@causes, collapse = ", ")
    ), call. = FALSE)
  list(cif = cif$cif, time_grid = cif$time_grid, model_key = cif$model_key)
}


#' Trapezoidal numerical integration
#'
#' @param x Numeric vector of x-values (must be sorted).
#' @param y Numeric vector of y-values (same length as `x`).
#'
#' @return A single numeric value.
#' @noRd
.trapezoidal_integration <- function(x, y) {
  if (length(x) != length(y))
    stop("`x` and `y` must have the same length.")
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}
