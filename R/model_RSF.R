#' @title Random Survival Forest for Competing Risks (RSF)
#' @description Registers the competing risks random survival forest via
#'   \code{randomForestSRC::rfsrc} into the \pkg{CompRisksBenchmaRk} registry.
#'
#' @details
#' Available under the key \code{"RSF"}.  Requires the **randomForestSRC**
#' package.  Hyperparameters (`ntree`, `mtry`, `nodesize`, `splitrule`) can
#' be tuned via the grid.
#'
#' @name model_RSF
#' @keywords internal
NULL

register_cr_model(
  key = "RSF",

  fit = function(data, formula, args = list()) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE))
      stop("Please install 'randomForestSRC'.")

    causes <- cr_causes(data, formula)
    fits   <- randomForestSRC::rfsrc(
      formula   = formula,
      data      = data,
      ntree     = cr_arg(args, "ntree",     500L),
      mtry      = cr_arg(args, "mtry",      NULL),
      nodesize  = cr_arg(args, "nodesize",  NULL),
      splitrule = cr_arg(args, "splitrule", "logrankCR")
    )
    structure(list(causes = causes, fits = fits),
              class = c("cr_model_rsf", "cr_model"))
  },

  predict_cif = function(fit_obj, newdata, times) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE))
      stop("Please install 'randomForestSRC'.")

    n       <- nrow(newdata)
    K       <- length(fit_obj$causes)
    out     <- array(0, dim = c(n, K, length(times)))
    pr      <- randomForestSRC::predict.rfsrc(fit_obj$fits, newdata = newdata)
    train_t <- pr$time.interest

    for (k in seq_len(K)) {
      cif_mat <- as.matrix(pr$cif[, , k, drop = FALSE][,, 1L])
      for (i in seq_len(n)) {
        out[i, k, ] <- stats::approx(
          x    = train_t,
          y    = cif_mat[i, ],
          xout = times,
          rule = 2,
          ties = "ordered"
        )$y
      }
    }
    out
  },

  info = function() list(
    name         = "Random Survival Forest (CR) \u2014 randomForestSRC::rfsrc",
    supports     = "CIF",
    needs_tuning = TRUE,
    default_grid = function() tibble::tibble(
      ntree     = 500L,
      nodesize  = list(NULL),
      mtry      = list(NULL),
      splitrule = "logrankCR"
    )
  )
)
