#' Compute the restricted mean lifetime (RMLT) from a CIF array
#'
#' Integrates each subject's cumulative incidence function (CIF) up to `maxT`
#' using trapezoidal integration, returning a per-subject RMLT for each
#' competing cause.  Predictions can be supplied either as a pre-computed CIF
#' list via `cif`, or as a fitted model object via `fit` (in which case the
#' CIF is computed internally using [predict_cif()]).  Exactly one of `cif` or
#' `fit` must be non-`NULL`.
#'
#' @param cr            A [cr_data()] object.
#' @param cif           A list as returned by [predict_cif()], with elements
#'   `$cif` (3-D numeric array `[n, K, Tm]`) and `$time_grid` (numeric vector).
#'   Mutually exclusive with `fit`.
#' @param fit           A fitted model object as returned by [fit_cr_model()].
#'   When supplied, the CIF is computed via [predict_cif()] before integration.
#'   Mutually exclusive with `cif`.
#' @param cif_time_grid Numeric vector of time points passed to [predict_cif()]
#'   when `fit` is non-`NULL`.  Must be provided when `fit` is supplied;
#'   ignored when `cif` is supplied directly (times are taken from
#'   `cif$time_grid`).
#' @param maxT          Upper integration limit. Defaults to the maximum
#'   observed event time in `cr`. A message is issued when the default is used.
#'
#' @return An `[n x K]` numeric matrix of per-subject RMLTs, with column names
#'   `"cause_1"`, `"cause_2"`, etc.
#' @export
compute_rmlt <- function(cr,
                         cif           = NULL,
                         fit           = NULL,
                         cif_time_grid = NULL,
                         maxT          = NULL) {
  
  if (!methods::is(cr, "cr_data"))
    stop("`cr` must be a cr_data object.", call. = FALSE)
  
  if (is.null(cif) && is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are non-NULL.",
         call. = FALSE)
  
  if (!is.null(maxT)) {
    if (!is.numeric(maxT) || length(maxT) != 1 || is.na(maxT) || maxT <= 0)
      stop("`maxT` must be a single positive numeric value.", call. = FALSE)
  }
  
  if (!is.null(fit)) {
    if (is.null(cif_time_grid))
      stop("`cif_time_grid` must be provided when `fit` is non-NULL.",
           call. = FALSE)
    cif <- predict_cif(fit, newdata = cr, time_grid = cif_time_grid)
  }
  
  # cif is a list(cif = [n,K,Tm] array, time_grid = numeric vector)
  time_grid <- cif$time_grid
  cif       <- cif$cif
  
  if (is.null(maxT)) {
    maxT <- max(cr@data[[cr@time_var]][cr@data[[cr@event_var]] != cr@cens_code])
    message(
      "`maxT` not supplied; defaulting to the maximum observed event time (",
      round(maxT, 4), ")."
    )
  }
  
  d     <- dim(cif)
  n     <- d[1]; K <- d[2]
  idx_t <- which(time_grid <= maxT)
  
  if (length(idx_t) < 2) {
    res <- matrix(0, nrow = n, ncol = K)
    colnames(res) <- paste0("cause_", cr@causes)
    return(res)
  }
  
  times_use <- time_grid[idx_t]
  res <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    cif_mat  <- matrix(cif[, k, idx_t, drop = FALSE], nrow = n)
    res[, k] <- apply(cif_mat, 1, .trapezoidal_integration, x = times_use)
  }
  colnames(res) <- paste0("cause_", cr@causes)
  res
}