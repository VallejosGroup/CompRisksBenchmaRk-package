#' @title Internal cross-validation utilities
#' @description Internal helpers for nested cross-validation workflows.
#' @name utils_cv
NULL


#' Pool a list of CIF arrays by taking the element-wise mean
#'
#' Used to combine CIF predictions across multiple imputed datasets.
#'
#' @param cif_list A list of 3-D numeric arrays, each of dimension
#'   `[n, K, Tm]`.
#'
#' @return A single 3-D array of the same dimensions.
#' @noRd
.pool_cifs_mean <- function(cif_list) {
  m      <- length(cif_list)
  pooled <- cif_list[[1]]
  if (m > 1)
    for (i in 2:m) pooled <- pooled + cif_list[[i]]
  pooled / m
}
