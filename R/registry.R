#' @title Competing Risks Model Registry
#'
#' @description
#' A lightweight registry that maps string keys to model objects.  Each
#' registered model must expose three functions:
#' \describe{
#'   \item{fit}{\code{fit(obj, args, \dots)} — trains the model on a \code{cr_data}
#'     object (see \code{cr_data()}). Grid parameters are passed via \code{args};
#'     any additional arguments for the underlying fitting function via \code{\dots}.}
#'   \item{predict_cif}{\code{predict_cif(fit_obj, newdata, time_grid)} — returns a
#'     3-D array of cumulative incidence function (CIF) predictions with
#'     dimensions \code{[n, K, |time_grid|]}.}
#'   \item{info}{\code{info()} — returns a named list with at least \code{name},
#'     \code{supports}, \code{needs_tuning}, and \code{default_grid}.}
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
#' A thin wrapper around a model's `predict_cif` function that extracts the
#' data frame from a [cr_data()] object (or accepts a plain data frame) and
#' coerces `time_grid` to numeric before dispatching.
#' The registered model is looked up from `fit_obj$model_key`.
#'
#' @param fit_obj A fitted model object returned by [fit_cr_model()].
#' @param newdata A [cr_data()] object or a plain data frame of new
#'   observations.
#' @param time_grid Numeric vector of evaluation times.
#'
#' @return A 3-D numeric array of dimensions `[n, K, length(time_grid)]`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_cr_model("FGR", cr)
#' cif <- predict_cif(fit, newdata = cr, time_grid = time_grid)
#' }
predict_cif <- function(fit_obj, newdata, time_grid) {
  mdl       <- get_cr_model(fit_obj$model_key)
  newdata   <- if (methods::is(newdata, "cr_data")) newdata@data
               else if (is.data.frame(newdata)) newdata
               else as.data.frame(newdata)
  time_grid <- as.numeric(time_grid)
  mdl$predict_cif(fit_obj, newdata, time_grid)
}


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


#' List all registered model keys
#'
#' @return A character vector of registered model keys (may be empty).
#' @export
#'
#' @examples
#' list_cr_models()
list_cr_models <- function() ls(envir = .cr_models)

