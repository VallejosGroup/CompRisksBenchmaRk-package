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

    x       <- as.data.frame(x)
    y_time  <- as.numeric(y_time)
    y_event <- as.integer(as.character(y_event))
    causes  <- sort(unique(y_event[y_event != 0L]))

    dat    <- data.frame(time = y_time, event = y_event, x,
                         check.names = FALSE)
    covars <- setdiff(names(dat), c("time", "event"))
    f_full <- stats::reformulate(covars, response = NULL)
    f_full <- stats::update(f_full, prodlim::Hist(time, event) ~ .)

    fits <- lapply(causes, function(k) {
      riskRegression::CSC(formula = f_full, data = dat, cause = k)
    })
    structure(list(causes = causes, fits = fits),
              class = c("cr_model_csCPH", "cr_model"))
  },

  predict_cif = function(fit_obj, x_new, times) {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Please install 'riskRegression'.")

    x_new <- as.data.frame(x_new)
    times <- as.numeric(times)
    n     <- nrow(x_new)
    K     <- length(fit_obj$causes)
    out   <- array(0, dim = c(n, K, length(times)))

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
