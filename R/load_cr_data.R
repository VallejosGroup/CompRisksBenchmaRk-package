#' Load a stored cr_data object from disk
#'
#' Reads a Parquet data file and a JSON metadata file that were written by
#' [cr_data()] with \code{store = TRUE} and reconstructs the original
#' \code{cr_data} S4 object.
#'
#' @param store_dir  Path to the directory containing the stored files.
#' @param store_name Optional character string prefix used when the files were
#'   saved.  If \code{NULL} (default), the function looks for
#'   \code{cr_data.parquet} and \code{cr_metadata.json}; otherwise it looks
#'   for \code{<store_name>_data.parquet} and
#'   \code{<store_name>_metadata.json}.
#'
#' @return A \code{cr_data} S4 object identical to the one that was stored.
#' @export
load_cr_data <- function(store_dir, store_name = NULL) {
  
  if (!is.character(store_dir) || length(store_dir) != 1 || !nzchar(store_dir))
    stop("`store_dir` must be a single non-empty string.", call. = FALSE)
  if (!dir.exists(store_dir))
    stop(sprintf("`store_dir` '%s' does not exist.", store_dir), call. = FALSE)
  
  prefix    <- if (!is.null(store_name)) store_name else "cr"
  data_path <- file.path(store_dir, paste0(prefix, "_data.parquet"))
  meta_path <- file.path(store_dir, paste0(prefix, "_metadata.json"))
  
  if (!file.exists(data_path))
    stop(sprintf("Data file not found: '%s'.", data_path), call. = FALSE)
  if (!file.exists(meta_path))
    stop(sprintf("Metadata file not found: '%s'.", meta_path), call. = FALSE)
  
  meta <- jsonlite::read_json(meta_path, simplifyVector = TRUE)
  data <- .apply_var_types(
    as.data.frame(arrow::read_parquet(data_path)),
    meta$var_types,
    meta$var_names
  )
  
  methods::new("cr_data",
               data         = data,
               causes       = as.integer(meta$causes),
               var_names    = meta$var_names,
               var_types    = meta$var_types,
               time_var     = meta$time_var,
               event_var    = meta$event_var,
               id_var       = meta$id_var,
               sort_by_time = isTRUE(meta$sort_by_time),
               time_offset  = as.numeric(meta$time_offset),
               cens_code    = as.integer(meta$cens_code)
  )
}