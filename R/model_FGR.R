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
      riskRegression::FGR(f_full, data = dat, cause = k)
    })
    structure(list(causes = causes, fits = fits),
              class = c("cr_model_fgr", "cr_model"))
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
      pj <- riskRegression::predictRisk(fit_obj$fits[[j]],
                                        newdata = x_new,
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
