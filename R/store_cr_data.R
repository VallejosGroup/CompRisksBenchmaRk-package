#' Store a cr_data object to disk
#'
#' Writes a \code{cr_data} object to disk as two files: a Parquet file
#' containing \code{cr@data} and a JSON file containing all metadata slots
#' (as returned by [cr_metadata()]).  The files can be read back with
#' [load_cr_data()].
#'
#' @param cr        A \code{cr_data} object.
#' @param store_dir Path to an existing directory where the files will be
#'   written.
#' @param store_name Optional character string used as a prefix for the output
#'   file names.  Files are named \code{<store_name>_data.parquet} and
#'   \code{<store_name>_metadata.json}.  If \code{NULL} (default), the prefix
#'   \code{"cr"} is used, giving \code{cr_data.parquet} and
#'   \code{cr_metadata.json}.
#' @param overwrite Logical; if \code{FALSE} (default), an error is raised
#'   when either output file already exists.  Set to \code{TRUE} to allow
#'   overwriting.
#'
#' @return Invisibly returns a named character vector with elements
#'   \code{data} and \code{metadata} giving the paths of the written files.
#' @export
store_cr_data <- function(cr,
                           store_dir,
                           store_name = NULL,
                           overwrite  = FALSE) {

  if (!methods::is(cr, "cr_data"))
    stop("`cr` must be a cr_data object.", call. = FALSE)
  if (!is.character(store_dir) || length(store_dir) != 1 || !nzchar(store_dir))
    stop("`store_dir` must be a single non-empty string.", call. = FALSE)
  if (!dir.exists(store_dir))
    stop(sprintf("`store_dir` '%s' does not exist.", store_dir), call. = FALSE)

  prefix    <- if (!is.null(store_name)) store_name else "cr"
  data_path <- file.path(store_dir, paste0(prefix, "_data.parquet"))
  meta_path <- file.path(store_dir, paste0(prefix, "_metadata.json"))

  existing <- c(data_path, meta_path)[file.exists(c(data_path, meta_path))]
  if (length(existing) > 0 && !overwrite)
    stop(
      "The following file(s) already exist and would be overwritten:\n",
      paste0("  ", existing, collapse = "\n"), "\n",
      "Re-run with `overwrite = TRUE` to allow overwriting.",
      call. = FALSE
    )

  arrow::write_parquet(cr_get_data(cr), data_path)
  jsonlite::write_json(cr_metadata(cr), path = meta_path,
                       pretty = TRUE, auto_unbox = TRUE)

  invisible(c(data = data_path, metadata = meta_path))
}
