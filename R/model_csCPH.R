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

  fit = function(data, time_var, event_var, args = list()) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    if (!requireNamespace("prodlim", quietly = TRUE))
      stop("Please install 'prodlim'.")

    causes  <- cr_causes(data, event_var)
    covars  <- setdiff(names(data), c(time_var, event_var))
    formula <- stats::reformulate(covars,
                 response = paste0("prodlim::Hist(", time_var, ", ", event_var, ")"))

    fits <- lapply(causes, function(k) {
      riskRegression::CSC(formula = formula, data = data, cause = k)
    })
    structure(list(causes = causes, fits = fits),
              class = c("cr_model_csCPH", "cr_model"))
  },

  predict_cif = function(fit_obj, newdata, time_grid) {
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
