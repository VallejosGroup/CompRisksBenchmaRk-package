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

  fit = function(data, time.var, event.var, args = list()) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")
    if (!requireNamespace("prodlim", quietly = TRUE))
      stop("Please install 'prodlim'.")

    causes  <- cr_causes(data, event.var)
    covars  <- setdiff(names(data), c(time.var, event.var))
    formula <- stats::reformulate(covars,
                 response = paste0("prodlim::Hist(", time.var, ", ", event.var, ")"))

    fits <- lapply(causes, function(k) {
      riskRegression::FGR(formula, data = data, cause = k)
    })
    structure(list(causes = causes, fits = fits),
              class = c("cr_model_fgr", "cr_model"))
  },

  predict_cif = function(fit_obj, newdata, times) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")

    n   <- nrow(newdata)
    K   <- length(fit_obj$causes)
    out <- array(0, dim = c(n, K, length(times)))
    for (j in seq_len(K)) {
      pj <- riskRegression::predictRisk(fit_obj$fits[[j]],
                                        newdata = newdata,
                                        times   = times)
      out[, j, ] <- as.matrix(pj)
    }
    out
  },

  info = function() list(
    name         = "Fine-Gray (riskRegression::FGR)",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble()
  )
)
