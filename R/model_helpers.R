#' @title Internal helpers shared across model registrations
#' @name model_helpers
#' @keywords internal
NULL


#' Extract competing causes from a data frame and formula
#'
#' Derives the sorted vector of competing causes (all non-zero event codes)
#' by extracting the event variable name from the formula response
#' (e.g. `prodlim::Hist(time, event)` or `survival::Surv(time, event)`)
#' and looking it up in `data`.
#'
#' @param data A data frame containing the event column.
#' @param formula A model formula whose LHS names the time and event variables.
#'
#' @return A sorted integer vector of cause codes.
#' @noRd
cr_causes <- function(data, formula) {
  lhs_vars  <- all.vars(formula[[2L]])   # e.g. c("time", "event")
  event_col <- lhs_vars[[2L]]            # second variable is always the event
  event     <- data[[event_col]]
  if (!is.integer(event))
    event <- as.integer(as.character(event))
  sort(unique(event[event != 0]))
}


#' Extract an argument from an `args` list with a default fallback
#'
#' @param args Named list of arguments (typically passed from the grid).
#' @param name Argument name (character scalar).
#' @param default Value to use when `name` is absent or `NULL` in `args`.
#'
#' @return `args[[name]]` if present and non-`NULL`, otherwise `default`.
#' @noRd
cr_arg <- function(args, name, default = NULL) {
  val <- args[[name]]
  if (is.null(val)) default else val
}
