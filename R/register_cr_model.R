#' Register a competing risks model
#'
#' Adds (or replaces) a model in the internal registry so it can be retrieved
#' by key and used inside [nested_cv_from_bench()].
#'
#' @param key A single non-empty character string used to identify the model.
#' @param fit A function with signature \code{fit(obj, args, \dots)} where \code{obj}
#'   is a \code{cr_data} object, \code{args} is a named list of grid-tuned
#'   arguments, and \code{\dots} accepts further arguments for the underlying
#'   fitting function.
#' @param predict_cif A function with signature
#'   \code{predict_cif(fit_obj, newdata, time_grid)} that returns a 3-D numeric
#'   array of dimensions \code{[n, K, Tm]}.
#' @param info A zero-argument function returning a named list with fields
#'   \code{name}, \code{supports}, \code{needs_tuning}, and \code{default_grid}.
#'
#' @return Invisibly `TRUE`.
#' @export
#'
#' @examples
#' register_cr_model(
#'   key = "my_model",
#'   fit = function(cr, args = list(), ...) {
#'     list(intercept = mean(cr@data[[cr@time_var]]), model_key = "my_model")
#'   },
#'   predict_cif = function(fit_obj, newdata, time_grid) {
#'     n <- nrow(newdata); K <- 1; Tm <- length(time_grid)
#'     array(0.1, dim = c(n, K, Tm))
#'   },
#'   info = function() list(name = "Trivial model", supports = "CIF",
#'                          needs_tuning = FALSE,
#'                          default_grid = function() tibble::tibble())
#' )
register_cr_model <- function(key, fit, predict_cif, info) {
  stopifnot(is.character(key), length(key) == 1, nzchar(key))
  stopifnot(is.function(fit), is.function(predict_cif), is.function(info))

  if (exists(key, envir = .cr_models, inherits = FALSE))
    warning("Overwriting registered model: ", key, call. = FALSE)

  assign(key, list(fit = fit, predict_cif = predict_cif, info = info),
         envir = .cr_models)

  invisible(TRUE)
}
