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
    
    inp    <- cr_prepare_inputs(x, time = y_time, event = y_event)
    x_cols <- names(inp$x)
    
    built  <- cr_build_formula(inp$x, inp$time, inp$event, "Hist")
    covars <- setdiff(names(built$dat), c("time", "event"))
    
    lambda <- cr_arg(args, "lambda_seq", NULL)
    if (is.null(lambda))
      stop("FGRP: lambda_seq must be provided via the grid.")
    
    fits <- lapply(inp$causes, function(k) {
      frm_k <- stats::as.formula(paste0(
        "fastcmprsk::Crisk(time, event, cencode=0, failcode=", k, ") ~ ",
        paste(covars, collapse = " + ")
      ))
      fp <- fastcmprsk::fastCrrp(
        formula     = frm_k,
        data        = built$dat,
        penalty     = cr_arg(args, "penalty",     "LASSO"),
        nlambda     = 1L,
        lambda      = as.numeric(lambda),
        standardize = cr_arg(args, "standardize", TRUE),
        alpha       = cr_arg(args, "alpha",        0)
      )
      if (!is.null(fp[["coef"]]) &&
          length(fp[["coef"]]) == length(x_cols)) {
        cf <- as.vector(fp[["coef"]])
        names(cf) <- x_cols
        fp[["coef"]] <- cf
      }
      list(fp = fp)
    })
    
    structure(list(causes = inp$causes, fits = fits, x_cols = x_cols),
              class = c("cr_model_fgrp", "cr_model"))
  },
  
  predict_cif = function(fit_obj, x_new, times) {
    p <- cr_init_cif_array(fit_obj, x_new, times)
    
    for (k in seq_len(p$K)) {
      fp    <- fit_obj$fits[[k]]$fp
      X     <- as.matrix(p$x_new[, !names(p$x_new) %in%
                                   c("time", "event", "row_id")])
      beta  <- as.vector(fp[["coef"]])
      bfitj <- fp[["breslowJump"]][[2L]]
      tt    <- as.vector(fp[["uftime"]])
      
      H0       <- cumsum(bfitj)
      H0_times <- stats::approx(x = tt, y = H0, xout = p$times,
                                method = "constant", f = 0, rule = 2)$y
      eta  <- drop(X %*% beta)
      lhat <- as.matrix(H0_times) %*% t(exp(eta))
      p$out[, k, ] <- t(1 - exp(-lhat))
    }
    p$out
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