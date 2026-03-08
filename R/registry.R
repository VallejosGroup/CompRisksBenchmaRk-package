#' @title Competing Risks Model Registry
#'
#' @description
#' A lightweight registry that maps string keys to model objects.  Each
#' registered model must expose three functions:
#' \describe{
#'   \item{fit}{`fit(data, time.var, event.var, args)` — trains the model on a
#'     data frame, using `time.var` and `event.var` to identify the time and
#'     event columns.}
#'   \item{predict_cif}{`predict_cif(fit_obj, newdata, times)` — returns a
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
#' @param fit A function with signature `fit(data, time.var, event.var, args)`
#'   where `data` is a data frame and `time.var` / `event.var` are character
#'   strings naming the respective columns.
#' @param predict_cif A function with signature
#'   `predict_cif(fit_obj, newdata, times)` that returns a 3-D numeric array
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
#'   fit = function(data, time.var, event.var, args = list()) {
#'     list(intercept = mean(data[[time.var]]))
#'   },
#'   predict_cif = function(fit_obj, newdata, times) {
#'     n <- nrow(newdata); K <- 1; Tm <- length(times)
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

  # Wrap fit so every returned object carries model_key automatically
  fit_wrapped <- function(data, time.var, event.var, args = list()) {
    fit_obj <- fit(data, time.var, event.var, args)
    fit_obj$model_key <- key
    fit_obj
  }

  assign(key, list(fit = fit_wrapped, predict_cif = predict_cif, info = info),
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
#' `newdata` to a data frame and `times` to numeric before dispatching.
#' The registered model is looked up from `fit_obj$model_key`, which is
#' stamped onto every fitted object automatically by [register_cr_model()].
#'
#' @param fit_obj A fitted model object returned by the model's `fit` function.
#' @param newdata Data frame of new observations.
#' @param times Numeric vector of evaluation times.
#'
#' @return A 3-D numeric array of dimensions `[n, K, length(times)]`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- get_cr_model("FGR")$fit(train, time.var = "time",
#'                                event.var = "event", args = list())
#' cif <- predict_cif(fit, newdata = test, times = times)
#' }
predict_cif <- function(fit_obj, newdata, times) {
  mdl     <- get_cr_model(fit_obj$model_key)
  newdata <- if (is.data.frame(newdata)) newdata else as.data.frame(newdata)
  times   <- as.numeric(times)
  mdl$predict_cif(fit_obj, newdata, times)
}


#' List all registered model keys
#'
#' @return A character vector of registered model keys (may be empty).
#' @export
#'
#' @examples
#' list_cr_models()
list_cr_models <- function() ls(envir = .cr_models)

