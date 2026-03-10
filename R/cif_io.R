#' @title CIF Array and Storage Utilities
#' @description Functions for pooling, saving, and storing CIF prediction
#'   arrays.
#' @name cif_io
NULL


#' Pool a list of CIF arrays by taking the element-wise mean
#'
#' Used to combine predictions across multiple imputed datasets.
#'
#' @param cif_list A list of 3-D numeric arrays, each of dimension
#'   `[n, K, Tm]`.
#'
#' @return A single 3-D array of the same dimensions.
#' @export
pool_cifs_mean <- function(cif_list) {
  m      <- length(cif_list)
  pooled <- cif_list[[1]]
  if (m > 1)
    for (i in 2:m) pooled <- pooled + cif_list[[i]]
  pooled / m
}



#' Build storage path for CIF outputs
#'
#' Creates the directory `out_dir/<model>/outer_<fold>/` and returns the path.
#'
#' @param out_dir  Root output directory.
#' @param model    Model key string.
#' @param fold     Outer fold index.
#'
#' @return The path (created on disk if needed).
#' @export
build_store_paths_r <- function(out_dir, model, fold) {
  p <- file.path(out_dir, model, sprintf("outer_%d", fold))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  p
}


#' Save CIF predictions to a Parquet file
#'
#' Flattens the 3-D CIF array into a long data frame with columns
#' `row_id`, `cause`, `time`, and `cif`, then writes it as a Parquet file.
#'
#' @param out_dir Directory in which `cif.parquet` is written.
#' @param cif     3-D array with dimensions `[n, K, Tm]`.
#' @param times   Numeric vector of evaluation times (length `Tm`).
#' @param row_ids Character vector of subject identifiers (length `n`).
#' @param causes  Integer vector of cause codes (length `K`).
#'
#' @return Invisibly `NULL`.
#' @export
save_cif_r <- function(out_dir, cif, times, row_ids, causes) {
  n  <- dim(cif)[1]; K <- dim(cif)[2]; Tm <- dim(cif)[3]
  df <- data.frame(
    row_id = rep(row_ids, each = K * Tm),
    cause  = rep(rep(causes, each = Tm), times = n),
    time   = rep(times, times = n * K),
    cif    = as.vector(cif)
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  arrow::write_parquet(df, file.path(out_dir, "cif.parquet"))
  invisible(NULL)
}