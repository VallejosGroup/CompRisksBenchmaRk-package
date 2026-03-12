#' @title Cause-Specific Cox Proportional Hazards Model (csCPH)
#' @description Registers the cause-specific Cox PH model via
#'   \code{riskRegression::CSC} into the \pkg{CompRisksBenchmaRk} registry.
#'
#' @details
#' Available under the key \code{"csCPH"}.  One CSC model is fitted per
#' competing cause; no hyperparameter tuning is required.
#'
#' @name model_csCPH
#' @keywords internal
NULL

register_cr_model(
  key = "csCPH",

  fit = function(obj, args = list(), ...) {
    if (!methods::is(obj, "cr_data"))
      stop("`obj` must be a cr_data object.", call. = FALSE)

    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    if (!requireNamespace("prodlim", quietly = TRUE))
      stop("Please install 'prodlim'.")

    formula <- stats::reformulate(
      obj@feature_vars,
      response = paste0("Hist(", obj@time_var, ", ", obj@event_var, ")")
    )
    environment(formula) <- asNamespace("prodlim")

    fits <- lapply(obj@causes, function(k) {
      riskRegression::CSC(formula = formula, data = obj@data, cause = k, ...)
    })
    list(causes = obj@causes, fits = fits, model_key = "csCPH")
  },

  predict_cif = function(fit_obj, newdata, time_grid) {
    if (!methods::is(newdata, "cr_data"))
      stop("`newdata` must be a cr_data object.", call. = FALSE)
    newdata <- newdata@data

    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")

    n   <- nrow(newdata)
    K   <- length(fit_obj$causes)
    out <- array(0, dim = c(n, K, length(time_grid)))
    for (j in seq_len(K)) {
      k  <- fit_obj$causes[j]
      pj <- riskRegression::predictRisk(
        fit_obj$fits[[j]],
        newdata       = newdata,
        times         = time_grid,
        cause         = k,
        product.limit = FALSE
      )
      out[, j, ] <- as.matrix(pj)
    }
    out
  },

  info = function() list(
    name         = "Cause-specific Cox PH \u2014 riskRegression::CSC",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble()
  )
)
