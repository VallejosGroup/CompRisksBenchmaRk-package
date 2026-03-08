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

  fit = function(data, time_var, event_var, args = list(), ...) {
    if (!requireNamespace("fastcmprsk", quietly = TRUE))
      stop("Please install 'fastcmprsk'.")

    causes <- cr_causes(data, event_var)
    x_cols <- setdiff(names(data), c(time_var, event_var, "row_id"))

    lambda <- args$lambda_seq
    if (is.null(lambda))
      stop("FGRP: lambda_seq must be provided via the grid.")

    fits <- lapply(causes, function(k) {
      frm_k <- stats::as.formula(paste0(
        "fastcmprsk::Crisk(", time_var, ", ", event_var,
        ", cencode=0, failcode=", k, ") ~ ",
        paste(x_cols, collapse = " + ")
      ))
      fp <- fastcmprsk::fastCrrp(
        formula     = frm_k,
        data        = data,
        penalty     = if (!is.null(args$penalty))     args$penalty     else "LASSO",
        nlambda     = 1L,
        lambda      = as.numeric(lambda),
        standardize = if (!is.null(args$standardize)) args$standardize else TRUE,
        alpha       = if (!is.null(args$alpha))       args$alpha       else 0,
        ...
      )
      if (!is.null(fp[["coef"]]) &&
          length(fp[["coef"]]) == length(x_cols)) {
        cf        <- as.vector(fp[["coef"]])
        names(cf) <- x_cols
        fp[["coef"]] <- cf
      }
      list(fp = fp)
    })

    structure(list(causes = causes, fits = fits, x_cols = x_cols, model_key = "FGRP"),
              class = c("cr_model_fgrp", "cr_model"))
  },

  predict_cif = function(fit_obj, newdata, time_grid) {
    n   <- nrow(newdata)
    K   <- length(fit_obj$causes)
    out <- array(0, dim = c(n, K, length(time_grid)))

    for (k in seq_len(K)) {
      fp    <- fit_obj$fits[[k]]$fp
      X     <- as.matrix(newdata[, fit_obj$x_cols, drop = FALSE])
      beta  <- as.vector(fp[["coef"]])
      bfitj <- fp[["breslowJump"]][[2L]]
      tt    <- as.vector(fp[["uftime"]])

      H0       <- cumsum(bfitj)
      H0_times <- stats::approx(x = tt, y = H0, xout = time_grid,
                                 method = "constant", f = 0, rule = 2)$y
      eta  <- drop(X %*% beta)
      lhat <- as.matrix(H0_times) %*% t(exp(eta))
      out[, k, ] <- t(1 - exp(-lhat))
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
