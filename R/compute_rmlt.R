#' Compute the restricted mean lifetime (RMLT) from a CIF array
#'
#' Integrates each subject's cumulative incidence function (CIF) up to `tau`
#' using trapezoidal integration, returning a per-subject RMLT for each
#' competing cause.  Predictions can be supplied either as a pre-computed CIF
#' array via `cif`, or as a fitted model object via `fit` (in which case the
#' CIF is computed internally using [predict_cif()]).  Exactly one of `cif` or
#' `fit` must be non-`NULL`.
#'
#' @param cr            A [cr_data()] object.
#' @param eval_times    Numeric vector of evaluation times (length `Tm`).
#' @param cif           3-D numeric array with dimensions `[n, K, Tm]`.
#'   Mutually exclusive with `fit`.
#' @param fit           A fitted model object as returned by [fit_cr_model()].
#'   When supplied, the CIF is computed via [predict_cif()] before integration.
#'   Mutually exclusive with `cif`.
#' @param cif_time_grid Numeric vector of time points passed to [predict_cif()]
#'   when `fit` is non-`NULL`.  Must be provided when `fit` is supplied;
#'   ignored when `cif` is supplied directly.
#' @param tau           Upper integration limit (default `max(eval_times)`).
#'
#' @return An `[n x K]` numeric matrix of per-subject RMLTs, with column names
#'   `"cause_1"`, `"cause_2"`, etc.
#' @export
compute_rmlt <- function(cr, eval_times,
                         cif           = NULL,
                         fit           = NULL,
                         cif_time_grid = NULL,
                         tau           = max(eval_times)) {

  if (!methods::is(cr, "cr_data"))
    stop("`cr` must be a cr_data object.", call. = FALSE)

  if (is.null(cif) && is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are non-NULL.",
         call. = FALSE)

  if (!is.null(fit)) {
    if (is.null(cif_time_grid))
      stop("`cif_time_grid` must be provided when `fit` is non-NULL.",
           call. = FALSE)
    cif <- predict_cif(fit, newdata = cr, time_grid = cif_time_grid)$cif
  }

  d     <- dim(cif)
  n     <- d[1]; K <- d[2]
  idx_t <- which(eval_times <= tau)

  if (length(idx_t) < 2) {
    res <- matrix(0, nrow = n, ncol = K)
    colnames(res) <- paste0("cause_", seq_len(K))
    return(res)
  }

  times_use <- eval_times[idx_t]
  res <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    cif_mat  <- matrix(cif[, k, idx_t, drop = FALSE], nrow = n)
    res[, k] <- apply(cif_mat, 1,
                      .trapezoidal_integration, x = times_use)
  }
  colnames(res) <- paste0("cause_", seq_len(K))
  res
}
