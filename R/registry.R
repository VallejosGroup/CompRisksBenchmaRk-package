#' @title Competing Risks Model Registry
#'
#' @description
#' A lightweight registry that maps string keys to model objects.  Each
#' registered model must expose three functions:
#' \describe{
#'   \item{fit}{`fit(x, y_time, y_event, args)` — trains the model and returns
#'     a model object.}
#'   \item{predict_cif}{`predict_cif(fit_obj, x_new, times)` — returns a
#'     3-D array of cumulative incidence function (CIF) predictions with
#'     dimensions `[n, K, |times|]`.}
#'   \item{info}{`info()` — returns a named list with at least `name`,
#'     `supports`, `needs_tuning`, and `default_grid`.}
#' }
#'
#' @name registry
NULL

# Internal environment that stores all registered models.
.cr_models <- new.env(parent = emptyenv())

#' Register a competing risks model
#'
#' Adds (or replaces) a model in the internal registry so it can be retrieved
#' by key and used inside [nested_cv_from_bench()].
#'
#' @param key A single non-empty character string used to identify the model.
#' @param fit A function with signature `fit(x, y_time, y_event, args)`.
#' @param predict_cif A function with signature
#'   `predict_cif(fit_obj, x_new, times)` that returns a 3-D numeric array
#'   of dimensions `[n, K, Tm]`.
#' @param info A zero-argument function returning a named list with fields
#'   `name`, `supports`, `needs_tuning`, and `default_grid`.
#'
#' @return Invisibly `TRUE`.
#' @export
#'
#' @examples
#' register_cr_model(
#'   key = "my_model",
#'   fit = function(x, y_time, y_event, args = list()) {
#'     list(intercept = mean(y_time))
#'   },
#'   predict_cif = function(fit_obj, x_new, times) {
#'     n <- nrow(x_new); K <- 1; Tm <- length(times)
#'     array(0.1, dim = c(n, K, Tm))
#'   },
#'   info = function() list(name = "Trivial model", supports = "CIF",
#'                          needs_tuning = FALSE,
#'                          default_grid = function() tibble::tibble())
#' )
register_cr_model <- function(key, fit, predict_cif, info) {
  stopifnot(is.character(key), length(key) == 1L, nzchar(key))
  stopifnot(is.function(fit), is.function(predict_cif), is.function(info))

  if (exists(key, envir = .cr_models, inherits = FALSE))
    warning("Overwriting registered model: ", key, call. = FALSE)

  assign(key, list(fit = fit, predict_cif = predict_cif, info = info),
         envir = .cr_models)

  invisible(TRUE)
}


#' Retrieve a registered competing risks model
#'
#' @param key Character string matching a previously registered model key.
#'
#' @return A named list with elements `fit`, `predict_cif`, and `info`.
#' @export
#'
#' @examples
#' \dontrun{
#' mdl <- get_cr_model("FGR")
#' mdl$info()
#' }
get_cr_model <- function(key) {
  m <- get0(key, envir = .cr_models, inherits = FALSE)
  if (is.null(m)) stop("Model not registered: ", key, call. = FALSE)
  m
}


#' Predict cumulative incidence functions for a registered model
#'
#' A thin wrapper around a model's `predict_cif` function that coerces
#' `x_new` to a data frame and `times` to numeric before dispatching,
#' so individual model implementations do not need to repeat this.
#'
#' @param key Character string identifying a registered model.
#' @param fit_obj A fitted model object returned by the model's `fit` function.
#' @param x_new Feature matrix or data frame for new observations.
#' @param times Numeric vector of evaluation times.
#'
#' @return A 3-D numeric array of dimensions `[n, K, length(times)]`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- get_cr_model("FGR")$fit(x_train, y_time, y_event)
#' cif <- predict_cif("FGR", fit, x_test, times)
#' }
predict_cif <- function(key, fit_obj, x_new, times) {
  mdl   <- get_cr_model(key)
  x_new <- if (is.data.frame(x_new)) x_new else as.data.frame(x_new)
  times <- as.numeric(times)
  mdl$predict_cif(fit_obj, x_new, times)
}


#' List all registered model keys
#'
#' @return A character vector of registered model keys (may be empty).
#' @export
#'
#' @examples
#' list_cr_models()
list_cr_models <- function() ls(envir = .cr_models)

