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

  fit = function(x, y_time, y_event, args = list()) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    if (!requireNamespace("prodlim", quietly = TRUE))
      stop("Please install 'prodlim'.")

    inp   <- cr_prepare_inputs(x, time = y_time, event = y_event)
    built <- cr_build_formula(inp$x, inp$time, inp$event, "Hist")

    fits <- lapply(inp$causes, function(k) {
      riskRegression::CSC(formula = built$formula, data = built$dat, cause = k)
    })
    structure(list(causes = inp$causes, fits = fits),
              class = c("cr_model_csCPH", "cr_model"))
  },

  predict_cif = function(fit_obj, x_new, times) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")

    n   <- nrow(x_new)
    K   <- length(fit_obj$causes)
    out <- array(0, dim = c(n, K, length(times)))
    for (j in seq_len(K)) {
      k  <- fit_obj$causes[j]
      pj <- riskRegression::predictRisk(
        fit_obj$fits[[j]],
        newdata       = x_new,
        times         = times,
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
