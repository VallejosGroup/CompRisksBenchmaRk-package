#' Predict cumulative incidence functions for a registered model
#'
#' A thin wrapper around a model's `predict_cif` function that validates
#' `newdata` as a [cr_data()] object and coerces `time_grid` to numeric
#' before dispatching.
#' The registered model is looked up from `fit_obj$model_key`.
#'
#' @param fit_obj A fitted model object returned by [fit_cr_model()].
#' @param newdata A [cr_data()] object of new observations.
#' @param time_grid Numeric vector of evaluation times.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`cif`}{A 3-D numeric array of dimensions `[n, K, length(time_grid)]`.}
#'     \item{`time_grid`}{The numeric vector of evaluation times used.}
#'     \item{`model_key`}{The character key of the model that produced the
#'       predictions, taken from `fit_obj$model_key`.}
#'     \item{`ids`}{The subject ID vector taken from
#'       `newdata@data[[newdata@id_var]]`, preserving the row order of
#'       `newdata`.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' fit    <- fit_cr_model("FGR", cr)
#' result <- predict_cif(fit, newdata = cr, time_grid = time_grid)
#' result$cif       # [n, K, Tm] array
#' result$time_grid # evaluation times
#' result$model_key # "FGR"
#' result$ids       # subject IDs aligned with result$cif rows
#' }
predict_cif <- function(fit_obj, newdata, time_grid) {
  if (!methods::is(newdata, "cr_data"))
    stop("`newdata` must be a cr_data object.", call. = FALSE)
  mdl       <- get_cr_model(fit_obj$model_key)
  time_grid <- as.numeric(time_grid)
  list(
    cif       = mdl$predict_cif(fit_obj, newdata, time_grid),
    time_grid = time_grid,
    model_key = fit_obj$model_key,
    ids       = newdata@data[[newdata@id_var]]
  )
}
