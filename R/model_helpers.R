#' @title Internal helpers shared across model registrations
#' @name model_helpers
#' @keywords internal
NULL


#' Extract competing causes from a data frame
#'
#' Derives the sorted vector of competing causes (all non-zero event codes)
#' from the named event column of a data frame, coercing to integer if needed.
#'
#' @param data A data frame containing the event column.
#' @param event.var Character string naming the event column in `data`.
#'
#' @return A sorted integer vector of cause codes.
#' @noRd
cr_causes <- function(data, event.var) {
  event <- data[[event.var]]
  if (!is.integer(event))
    event <- as.integer(as.character(event))
  sort(unique(event[event != 0]))
}
