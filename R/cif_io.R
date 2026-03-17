#' @title CIF Array and Storage Utilities
#' @description Functions for saving and storing CIF prediction arrays.
#' @name cif_io
NULL



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
#' @param out_dir       Directory in which `cif.parquet` is written.
#' @param cif           3-D array with dimensions `[n, K, Tm]`.
#' @param cif_time_grid Numeric vector of evaluation times (length `Tm`).
#' @param row_ids       Character vector of subject identifiers (length `n`).
#' @param causes        Integer vector of cause codes (length `K`).
#'
#' @return Invisibly `NULL`.
#' @export
save_cif_r <- function(out_dir, cif, cif_time_grid, row_ids, causes) {
  n  <- dim(cif)[1]; K <- dim(cif)[2]; Tm <- dim(cif)[3]
  df <- data.frame(
    row_id = rep(row_ids, each = K * Tm),
    cause  = rep(rep(causes, each = Tm), times = n),
    time   = rep(cif_time_grid, times = n * K),
    cif    = as.vector(cif)
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  arrow::write_parquet(df, file.path(out_dir, "cif.parquet"))
  invisible(NULL)
}