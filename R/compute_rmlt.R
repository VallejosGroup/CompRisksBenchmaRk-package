#' Compute the restricted mean lifetime (RMLT) from a CIF array
#'
#' Integrates each subject's CIF up to `tau` using trapezoidal integration.
#'
#' @param out   3-D CIF array with dimensions `[n, K, Tm]`.
#' @param times Numeric vector of evaluation times (length `Tm`).
#' @param tau   Upper integration limit (default `max(times)`).
#'
#' @return An `[n x K]` numeric matrix of per-subject RMLTs.
#' @export
compute_rmlt <- function(out, times, tau = max(times)) {
  d     <- dim(out)
  n     <- d[1]; K <- d[2]
  idx_t <- which(times <= tau)

  if (length(idx_t) < 2) {
    res <- matrix(0, nrow = n, ncol = K)
    colnames(res) <- paste0("cause_", seq_len(K))
    return(res)
  }

  times_use <- times[idx_t]
  res <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    cif_mat  <- matrix(out[, k, idx_t, drop = FALSE], nrow = n)
    res[, k] <- apply(cif_mat, 1,
                      trapezoidal.integration, x = times_use)
  }
  colnames(res) <- paste0("cause_", seq_len(K))
  res
}
