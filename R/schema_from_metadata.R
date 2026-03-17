#' Convert cr_metadata to legacy schema format
#'
#' Converts the metadata list returned by [cr_metadata()] into the schema
#' format previously produced by `infer_data_types()`, for backwards
#' compatibility with pipelines that read `data_types.json`.
#'
#' @param metadata A named list as returned by [cr_metadata()].
#'
#' @return A named list with elements `version`, `core_cols`,
#'   `feature_cols`, and `types`, compatible with the legacy
#'   `data_types.json` format.
#' @export
schema_from_metadata <- function(metadata) {
  core_cols     <- c(metadata$time_var, metadata$event_var, metadata$id_var)
  feature_vars  <- metadata$feature_vars
  feature_types <- metadata$feature_types

  types <- as.list(stats::setNames(
    c(rep(NA_character_, length(core_cols)), feature_types),
    c(core_cols, feature_vars)
  ))

  list(
    version      = 1,
    core_cols    = core_cols,
    feature_cols = feature_vars,
    types        = types
  )
}
