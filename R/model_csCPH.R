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
    
    inp   <- cr_prepare_inputs(x, y_time, y_event)
    built <- cr_build_formula(inp$x, inp$y_time, inp$y_event, "Hist")
    
    fits <- lapply(inp$causes, function(k) {
      riskRegression::CSC(formula = built$formula, data = built$dat, cause = k)
    })
    structure(list(causes = inp$causes, fits = fits),
              class = c("cr_model_csCPH", "cr_model"))
  },
  
  predict_cif = function(fit_obj, x_new, times) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    
    p <- cr_init_cif_array(fit_obj, x_new, times)
    for (j in seq_len(p$K)) {
      k  <- fit_obj$causes[j]
      pj <- riskRegression::predictRisk(
        fit_obj$fits[[j]],
        newdata       = p$x_new,
        times         = p$times,
        cause         = k,
        product.limit = FALSE
      )
      p$out[, j, ] <- as.matrix(pj)
    }
    p$out
  },
  
  info = function() list(
    name         = "Cause-specific Cox PH \u2014 riskRegression::CSC",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble()
  )
)