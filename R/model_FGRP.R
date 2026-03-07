#' @title Penalised Fine-Gray Regression Model (FGRP)
#' @description Registers the LASSO-penalised Fine-Gray model via
#'   \code{fastcmprsk::fastCrrp} into the \pkg{CompRisksBenchmaRk} registry.
#'
#' @details
#' Available under the key \code{"FGRP"}.  Requires the **fastcmprsk**
#' package.  This model requires hyperparameter tuning (`needs_tuning = TRUE`);
#' the default tuning grid sweeps 50 log-spaced \eqn{\lambda} values from 1
#' down to \eqn{10^{-4}}.
#'
#' @name model_FGRP
#' @keywords internal
NULL

register_cr_model(
  key = "FGRP",

  fit = function(x, y_time, y_event, args = list()) {
    if (!requireNamespace("fastcmprsk", quietly = TRUE))
      stop("Please install 'fastcmprsk'.")

    x       <- as.data.frame(x)
    x_cols  <- names(x)
    y_time  <- as.numeric(y_time)
    y_event <- as.integer(as.character(y_event))
    causes  <- sort(unique(y_event[y_event != 0L]))

    dat     <- data.frame(time = y_time, event = y_event, x,
                          check.names = FALSE)
    covars  <- setdiff(names(dat), c("time", "event"))

    penalty     <- if (!is.null(args$penalty))     args$penalty     else "LASSO"
    standardize <- if (!is.null(args$standardize)) args$standardize else TRUE
    alpha       <- if (!is.null(args$alpha))       args$alpha       else 0
    lambda      <- if (!is.null(args$lambda_seq))
      as.numeric(args$lambda_seq)
    else
      stop("FGRP: lambda_seq must be provided via the grid.")

    fits <- lapply(causes, function(k) {
      frm_k <- stats::as.formula(paste0(
        "fastcmprsk::Crisk(time, event, cencode=0, failcode=", k, ") ~ ",
        paste(covars, collapse = " + ")
      ))
      fp <- fastcmprsk::fastCrrp(
        formula     = frm_k,
        data        = dat,
        penalty     = penalty,
        nlambda     = 1L,
        lambda      = lambda,
        standardize = standardize,
        alpha       = alpha
      )
      if (!is.null(fp[["coef"]]) &&
          length(fp[["coef"]]) == length(x_cols)) {
        cf <- as.vector(fp[["coef"]])
        names(cf) <- x_cols
        fp[["coef"]] <- cf
      }
      list(fp = fp)
    })

    structure(list(causes = causes, fits = fits, x_cols = x_cols),
              class = c("cr_model_fgrp", "cr_model"))
  },

  predict_cif = function(fit_obj, x_new, times) {
    x_new <- as.data.frame(x_new)
    times <- as.numeric(times)
    n     <- nrow(x_new)
    K     <- length(fit_obj$causes)
    out   <- array(0, dim = c(n, K, length(times)))

    for (k in seq_len(K)) {
      fp    <- fit_obj$fits[[k]]$fp
      X     <- as.matrix(x_new[, !names(x_new) %in%
                                 c("time", "event", "row_id")])
      beta  <- as.vector(fp[["coef"]])
      bfitj <- fp[["breslowJump"]][[2L]]
      tt    <- as.vector(fp[["uftime"]])

      H0       <- cumsum(bfitj)
      H0_times <- stats::approx(x = tt, y = H0, xout = times,
                                 method = "constant", f = 0, rule = 2)$y
      eta  <- drop(X %*% beta)
      lhat <- as.matrix(H0_times) %*% t(exp(eta))
      cif  <- t(1 - exp(-lhat))
      out[, k, ] <- cif
    }
    out
  },

  info = function() list(
    name         = "Penalized Fine-Gray \u2014 fastcmprsk::fastCrrp",
    supports     = "CIF",
    needs_tuning = TRUE,
    default_grid = function() tibble::tibble(
      penalty    = "LASSO",
      lambda_seq = exp(seq(log(1), log(1e-4), length.out = 50))
    )
  )
)
