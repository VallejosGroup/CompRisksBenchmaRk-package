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

  fit = function(x, y_time, y_event, args = list()) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE))
      stop("Please install 'randomForestSRC'.")

    x       <- as.data.frame(x)
    y_time  <- as.numeric(y_time)
    y_event <- as.integer(as.character(y_event))
    causes  <- sort(unique(y_event[y_event != 0L]))

    ntree     <- if (!is.null(args$ntree))     args$ntree     else 500L
    mtry      <- if (!is.null(args$mtry))      args$mtry      else NULL
    nodesize  <- if (!is.null(args$nodesize))  args$nodesize  else NULL
    splitrule <- if (!is.null(args$splitrule)) args$splitrule else "logrankCR"

    dat    <- data.frame(time = y_time, event = y_event, x,
                         check.names = FALSE)
    covars <- setdiff(names(dat), c("time", "event"))
    f_full <- stats::reformulate(covars, response = NULL)
    f_full <- stats::update(f_full, survival::Surv(time, event) ~ .)

    fits <- randomForestSRC::rfsrc(
      formula   = f_full,
      data      = dat,
      ntree     = ntree,
      mtry      = mtry,
      nodesize  = nodesize,
      splitrule = splitrule
    )
    structure(list(causes = causes, fits = fits),
              class = c("cr_model_rsf", "cr_model"))
  },

  predict_cif = function(fit_obj, x_new, times) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE))
      stop("Please install 'randomForestSRC'.")

    x_new <- as.data.frame(x_new)
    times <- as.numeric(times)
    pr    <- randomForestSRC::predict.rfsrc(fit_obj$fit,
                                            newdata = x_new)
    train_t <- pr$time.interest
    n       <- nrow(x_new)
    K       <- length(fit_obj$causes)
    out     <- array(0, dim = c(n, K, length(times)))

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
