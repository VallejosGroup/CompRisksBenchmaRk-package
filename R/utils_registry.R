#' @title Internal registry utilities
#' @description Internal environment storing registered models, and helper
#'   functions shared across model registrations.
#' @name utils_registry
#' @keywords internal
NULL

# Internal environment that stores all registered models.
.cr_models <- new.env(parent = emptyenv())


#' Extract competing causes from a data frame
#'
#' Derives the sorted vector of competing causes (all non-zero event codes)
#' from the named event column of a data frame, coercing to integer if needed.
#'
#' @param data A data frame containing the event column.
#' @param event_var Character string naming the event column in `data`.
#'
#' @return A sorted integer vector of cause codes.
#' @noRd
cr_causes <- function(data, event_var) {
  event <- data[[event_var]]
  if (!is.integer(event))
    event <- as.integer(as.character(event))
  sort(unique(event[event != 0]))
}
