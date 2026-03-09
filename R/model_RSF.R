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

  fit = function(obj, args = list(), ...) {
    if (!methods::is(obj, "cr_data"))
      stop("`obj` must be a cr_data object.", call. = FALSE)

    if (!requireNamespace("randomForestSRC", quietly = TRUE))
      stop("Please install 'randomForestSRC'.")

    formula <- stats::reformulate(
      obj@covars$covars_names,
      response = paste0("survival::Surv(", obj@time_var, ", ", obj@event_var, ")")
    )

    fits <- randomForestSRC::rfsrc(
      formula   = formula,
      data      = obj@data,
      ntree     = if (!is.null(args$ntree))     args$ntree     else 500,
      mtry      = if (!is.null(args$mtry))      args$mtry      else NULL,
      nodesize  = if (!is.null(args$nodesize))  args$nodesize  else NULL,
      splitrule = if (!is.null(args$splitrule)) args$splitrule else "logrankCR",
      ...
    )
    list(causes = obj@causes, fits = fits, model_key = "RSF")
  },

  predict_cif = function(fit_obj, newdata, time_grid) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE))
      stop("Please install 'randomForestSRC'.")

    n       <- nrow(newdata)
    K       <- length(fit_obj$causes)
    out     <- array(0, dim = c(n, K, length(time_grid)))
    pr      <- randomForestSRC::predict.rfsrc(fit_obj$fits, newdata = newdata)
    train_t <- pr$time.interest

    for (k in seq_len(K)) {
      cif_mat <- as.matrix(pr$cif[, , k, drop = FALSE][,, 1])
      for (i in seq_len(n)) {
        out[i, k, ] <- stats::approx(
          x    = train_t,
          y    = cif_mat[i, ],
          xout = time_grid,
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
      ntree     = 500,
      nodesize  = list(NULL),
      mtry      = list(NULL),
      splitrule = "logrankCR"
    )
  )
)
