#' @title Competing Risks Model Registry
#'
#' @description
#' A lightweight registry that maps string keys to model objects.  Each
#' registered model must expose three functions:
#' \describe{
#'   \item{fit}{`fit(data, time_var, event_var, args, ...)` — trains the model on a
#'     data frame, using `time_var` and `event_var` to identify the time and
#'     event columns. Grid parameters are passed via `args`; any additional
#'     arguments for the underlying fitting function can be passed via `...`.}
#'   \item{predict_cif}{`predict_cif(fit_obj, newdata, time_grid)` — returns a
#'     3-D array of cumulative incidence function (CIF) predictions with
#'     dimensions `[n, K, |time_grid|]`.}
#'   \item{info}{`info()` — returns a named list with at least `name`,
#'     `supports`, `needs_tuning`, and `default_grid`.}
#' }
#'
#' @name registry
NULL

#' Register a competing risks model
#'
#' Adds (or replaces) a model in the internal registry so it can be retrieved
#' by key and used inside [nested_cv_from_bench()].
#'
#' @param key A single non-empty character string used to identify the model.
#' @param fit A function with signature `fit(data, time_var, event_var, args, ...)`
#'   where `data` is a data frame, `time_var` / `event_var` are character
#'   strings naming the respective columns, `args` is a named list of
#'   grid-tuned arguments, and `...` accepts any further arguments passed
#'   directly to the underlying fitting function.
#' @param predict_cif A function with signature
#'   `predict_cif(fit_obj, newdata, time_grid)` that returns a 3-D numeric array
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
#'   fit = function(data, time_var, event_var, args = list(), ...) {
#'     list(intercept = mean(data[[time_var]]), model_key = "my_model")
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
#' `newdata` to a data frame and `time_grid` to numeric before dispatching.
#' The registered model is looked up from `fit_obj$model_key`, which each
#' model's `fit` function stores in the returned object via `structure()`.
#'
#' @param fit_obj A fitted model object returned by the model's `fit` function.
#' @param newdata Data frame of new observations.
#' @param time_grid Numeric vector of evaluation time_grid.
#'
#' @return A 3-D numeric array of dimensions `[n, K, length(time_grid)]`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- get_cr_model("FGR")$fit(train, time_var = "time",
#'                                event_var = "event", args = list())
#' cif <- predict_cif(fit, newdata = test, time_grid = time_grid)
#' }
predict_cif <- function(fit_obj, newdata, time_grid) {
  mdl     <- get_cr_model(fit_obj$model_key)
  newdata <- if (is.data.frame(newdata)) newdata else as.data.frame(newdata)
  time_grid   <- as.numeric(time_grid)
  mdl$predict_cif(fit_obj, newdata, time_grid)
}


#' List all registered model keys
#'
#' @return A character vector of registered model keys (may be empty).
#' @export
#'
#' @examples
#' list_cr_models()
list_cr_models <- function() ls(envir = .cr_models)

