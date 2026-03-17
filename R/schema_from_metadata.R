#' Convert cr_metadata to legacy schema format
#'
#' Converts the metadata list returned by [cr_metadata()] into the legacy
#' schema format, for backwards compatibility with pipelines that read
#' `data_types.json`.
#'
#' @param metadata A named list as returned by [cr_metadata()].
#'
#' @return A named list with elements `version`, `core_cols`,
#'   `feature_cols`, and `types`, compatible with the legacy
#'   `data_types.json` format.
#' @export
schema_from_metadata <- function(metadata) {
  core_cols    <- c(metadata$time_var, metadata$event_var, metadata$id_var)
  feature_vars <- setdiff(metadata$var_names, c(metadata$time_var, metadata$event_var, metadata$id_var))

  # var_types is now a named vector covering core + feature cols
  types <- as.list(metadata$var_types)

  list(
    version      = 1,
    core_cols    = core_cols,
    feature_cols = feature_vars,
    types        = types
  )
}
