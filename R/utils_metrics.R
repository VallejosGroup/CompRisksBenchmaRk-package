#' @title Internal metric utilities
#' @description Internal helper functions shared across metric computations.
#' @name utils_metrics
#' @keywords internal
NULL


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
