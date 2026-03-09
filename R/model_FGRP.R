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
#' Factor covariates are expanded to dummy columns via \code{model.matrix}
#' before fitting; the terms reference is stored on the fit object for aligned
#' expansion at predict time.
#'
#' @name model_FGRP
#' @keywords internal
NULL

register_cr_model(
  key = "FGRP",

  fit = function(obj, args = list(), ...) {
    if (!methods::is(obj, "cr_data"))
      stop("`obj` must be a cr_data object.", call. = FALSE)

    if (!requireNamespace("fastcmprsk", quietly = TRUE))
      stop("Please install 'fastcmprsk'.")

    lambda <- args$lambda_seq
    if (is.null(lambda))
      stop("FGRP: lambda_seq must be provided via the grid.")

    X         <- .make_model_matrix(obj)
    terms_ref <- list(tt = attr(X, "tt"), col_names = attr(X, "col_names"))
    x_cols    <- terms_ref$col_names

    # fastCrrp requires a plain data frame with the expanded columns
    fit_data <- as.data.frame(X)
    fit_data[[obj@time_var]]  <- obj@data[[obj@time_var]]
    fit_data[[obj@event_var]] <- obj@data[[obj@event_var]]

    fits <- lapply(obj@causes, function(k) {
      frm_k <- stats::as.formula(paste0(
        "fastcmprsk::Crisk(", obj@time_var, ", ", obj@event_var,
        ", cencode=0, failcode=", k, ") ~ ",
        paste(x_cols, collapse = " + ")
      ))
      fp <- fastcmprsk::fastCrrp(
        formula     = frm_k,
        data        = fit_data,
        penalty     = if (!is.null(args$penalty))     args$penalty     else "LASSO",
        nlambda     = 1,
        lambda      = as.numeric(lambda),
        standardize = if (!is.null(args$standardize)) args$standardize else TRUE,
        alpha       = if (!is.null(args$alpha))       args$alpha       else 0,
        ...
      )
      cf        <- as.vector(fp[["coef"]])
      names(cf) <- x_cols
      fp[["coef"]] <- cf
      list(fp = fp)
    })

    list(causes = obj@causes, fits = fits, x_cols = x_cols,
         terms_ref = terms_ref, model_key = "FGRP")
  },

  predict_cif = function(fit_obj, newdata, time_grid) {
    if (!methods::is(newdata, "cr_data"))
      stop("`newdata` must be a cr_data object.", call. = FALSE)

    X   <- .make_model_matrix(newdata, fit_obj)
    n   <- nrow(X)
    K   <- length(fit_obj$causes)
    out <- array(0, dim = c(n, K, length(time_grid)))

    for (k in seq_len(K)) {
      fp    <- fit_obj$fits[[k]]$fp
      beta  <- as.vector(fp[["coef"]])
      bfitj <- fp[["breslowJump"]][[2]]
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
