#' @title Fine-Gray Regression Model (FGR)
#' @description Registers the Fine-Gray subdistribution hazard model via
#'   \code{riskRegression::FGR} into the \pkg{CompRisksBenchmaRk} registry.
#'   One model is fitted per competing cause.
#'
#' @details
#' This registration runs automatically when the package is loaded.
#' The model is available under the key \code{"FGR"}.
#'
#' **Required packages:** `riskRegression`, `prodlim`.
#'
#' **Hyperparameters:** None (`needs_tuning = FALSE`).
#'
#' @name model_FGR
#' @keywords internal
NULL

register_cr_model(
  key = "FGR",
  
  fit = function(x, y_time, y_event, args = list()) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    if (!requireNamespace("prodlim", quietly = TRUE))
      stop("Please install 'prodlim'.")
    
    inp    <- cr_prepare_inputs(x, time = y_time, event = y_event)
    built  <- cr_build_formula(inp$x, inp$time, inp$event, "Hist")
    
    fits <- lapply(inp$causes, function(k) {
      riskRegression::FGR(built$formula, data = built$dat, cause = k)
    })
    structure(list(causes = inp$causes, fits = fits),
              class = c("cr_model_fgr", "cr_model"))
  },
  
  predict_cif = function(fit_obj, x_new, times) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    
    p <- cr_init_cif_array(fit_obj, x_new, times)
    for (j in seq_len(p$K)) {
      pj <- riskRegression::predictRisk(fit_obj$fits[[j]],
                                        newdata = p$x_new,
                                        times   = p$times)
      p$out[, j, ] <- as.matrix(pj)
    }
    p$out
  },
  
  info = function() list(
    name         = "Fine-Gray (riskRegression::FGR)",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble()
  )
)