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
#'   must be `NULL` when `cif` is supplied directly.
#' @param maxT          Upper integration limit. Defaults to the maximum
#'   observed event time in `cr`. A message is issued when the default is used.
#'
#' @return An `[n x K]` numeric matrix of per-subject RMLTs, with column names
#'   `"cause_<k>"` where `<k>` is the actual cause code from `cr@causes`.
#' @export
compute_rmlt <- function(cr,
                         cif           = NULL,
                         fit           = NULL,
                         cif_time_grid = NULL,
                         maxT          = NULL) {
  
  .check_cr(cr)
  
  cif      <- .resolve_cif(cif, fit, cif_time_grid, cr)
  unpacked <- .validate_and_unpack_cif(cif, cr)
  time_grid <- unpacked$time_grid
  cif       <- unpacked$cif
  
  if (!is.null(maxT)) {
    if (!is.numeric(maxT) || length(maxT) != 1 || is.na(maxT) || maxT <= 0)
      stop("`maxT` must be a single positive numeric value.", call. = FALSE)
  }
  
  if (is.null(maxT)) {
    maxT <- max(cr@data[[cr@time_var]][cr@data[[cr@event_var]] != cr@cens_code])
    message(
      "`maxT` not supplied; defaulting to the maximum observed event time (",
      round(maxT, 4), ")."
    )
  }
  
  n     <- nrow(cif)
  K     <- dim(cif)[2]
  idx_t <- which(time_grid <= maxT)
  
  if (length(idx_t) < 2) {
    res <- matrix(0, nrow = n, ncol = K)
    colnames(res) <- paste0("cause_", cr@causes)
    return(res)
  }
  
  times_use <- time_grid[idx_t]
  res <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    cif_mat  <- cif[, k, idx_t, drop = TRUE]
    res[, k] <- apply(cif_mat, 1, .trapezoidal_integration, x = times_use)
  }
  colnames(res) <- paste0("cause_", cr@causes)
  res
}