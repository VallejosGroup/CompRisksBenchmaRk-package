#' Fit a registered competing risks model
#'
#' A convenience wrapper around [get_cr_model()] that retrieves the model by
#' key and calls its \code{fit} function directly, without needing to call
#' \code{get_cr_model(key)$fit(...)}.
#'
#' @param model_key A single non-empty character string identifying a
#'   registered model (see [list_cr_models()]).
#' @param obj A \code{cr_data} object produced by [cr_data()].
#' @param args A named list of hyperparameter arguments passed to the model's
#'   \code{fit} function.
#' @param ... Additional arguments forwarded to the underlying fitting
#'   function.
#'
#' @return A fitted model object (a named list containing at least
#'   \code{model_key} and \code{causes}).
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_cr_model("FGR", obj)
#' }
fit_cr_model <- function(model_key, obj, args = list(), ...) {
  get_cr_model(model_key)$fit(obj = obj, args = args, ...)
}
