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


#' Build or align a model matrix for FGRP
#'
#' At fit time (\code{fit_obj = NULL}): expands all covariates to a numeric
#' matrix via \code{model.matrix}, drops zero-variance columns, and stores
#' the \code{terms} object and final column names as attributes for later
#' alignment.
#'
#' At predict time (\code{fit_obj} supplied): re-expands \code{newdata} using
#' the stored \code{terms} reference, aligns columns to those seen at fit time
#' (filling missing columns with 0), and returns the numeric matrix.
#'
#' @param obj     A \code{cr_data} object.
#' @param fit_obj At predict time, the fitted FGRP list object (must contain
#'   \code{terms_ref} with elements \code{tt} and \code{col_names}).
#'   Pass \code{NULL} at fit time.
#'
#' @return A numeric matrix.  At fit time the matrix carries two attributes:
#'   \code{tt} (the \code{terms} object) and \code{col_names} (character
#'   vector of retained column names).
#' @noRd
.make_model_matrix <- function(obj, fit_obj = NULL) {
  data <- obj@data

  if (is.null(fit_obj)) {
    # --- fit time ---
    rhs     <- paste(obj@covars$covars_names, collapse = " + ")
    formula <- stats::as.formula(paste("~", rhs))
    mf      <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
    tt      <- attr(mf, "terms")
    X       <- stats::model.matrix(tt, data = mf)[, -1, drop = FALSE]  # drop intercept

    # drop zero-variance columns
    keep <- apply(X, 2, stats::var) > 0
    X    <- X[, keep, drop = FALSE]

    attr(X, "tt")        <- tt
    attr(X, "col_names") <- colnames(X)
    X

  } else {
    # --- predict time ---
    tt        <- fit_obj$terms_ref$tt
    col_names <- fit_obj$terms_ref$col_names

    mf <- stats::model.frame(tt, data = data,
                             na.action    = stats::na.pass,
                             xlev         = attr(tt, "xlev"))
    X  <- stats::model.matrix(tt, data = mf)[, -1, drop = FALSE]

    # align to fit-time columns, filling any missing with 0
    missing_cols <- setdiff(col_names, colnames(X))
    extra        <- matrix(0, nrow = nrow(X),
                           ncol = length(missing_cols),
                           dimnames = list(NULL, missing_cols))
    X <- cbind(X, extra)[, col_names, drop = FALSE]
    X
  }
}
