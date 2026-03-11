#' Calibration plots for competing risks
#'
#' Produces pseudo-value calibration plots using LOESS or nearest-neighbour
#' smoothing for all competing causes.  Returns plot objects, smoothed
#' calibration frames, and numerical calibration measures in a format
#' consistent with [compute_metrics()].
#'
#' @param cr            A [cr_data()] object providing outcomes and metadata.
#' @param cif           A list as returned by [predict_cif()], with elements
#'   `$cif` (3-D numeric array `[n, K, Tm]`), `$time_grid` (numeric vector),
#'   and optionally `$model_key` (used for plot titles).
#'   Mutually exclusive with `fit`.
#' @param fit           A fitted model object as returned by [fit_cr_model()].
#'   When supplied, the CIF is computed via [predict_cif()] before plotting.
#'   Mutually exclusive with `cif`.
#' @param cif_time_grid Numeric vector of time points passed to [predict_cif()]
#'   when `fit` is non-`NULL`. Must be `NULL` when `cif` is supplied directly.
#' @param tau           Numeric scalar or vector of evaluation times.
#' @param loess_smoothing If `TRUE` (default), use LOESS smoothing;
#'   otherwise use `prodlim::meanNeighbors()`.
#' @param bandwidth     Optional bandwidth for nearest-neighbour smoothing.
#' @param graph         Logical; produce ggplot objects? (default `TRUE`).
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{`graphs`}{Named list (one entry per cause, e.g. `"cause_1"`) of
#'       plot lists, each with one ggplot per element of `tau`.  `NULL` when
#'       `graph = FALSE`.}
#'     \item{`values`}{Named list (one entry per cause) of smoothed calibration
#'       data frames, each with one data frame per element of `tau`.}
#'     \item{`calib_measures`}{Data frame with columns `cause`, `tau`, `ICI`,
#'       `E50`, `E90`, `Emax`, `RSB`.  One row per cause per evaluation time.}
#'     \item{`OE_summary`}{Data frame with columns `cause`, `OE`, `lower`,
#'       `upper` (95\% CI).  One row per cause, evaluated at the last element
#'       of `tau`.}
#'   }
#' @export
CalibrationPlot <- function(cr,
                             cif           = NULL,
                             fit           = NULL,
                             cif_time_grid = NULL,
                             tau,
                             loess_smoothing = TRUE,
                             bandwidth       = NULL,
                             graph           = TRUE) {
  .check_cr(cr)

  if (missing(tau) || !is.numeric(tau) || length(tau) == 0)
    stop("`tau` must be a non-empty numeric vector.", call. = FALSE)
  if (anyNA(tau) || is.unsorted(tau))
    stop("`tau` must be sorted and contain no NA values.", call. = FALSE)

  cif_in   <- .resolve_cif(cif, fit, cif_time_grid, cr)
  unpacked <- .validate_and_unpack_cif(cif_in, cr)
  cif_arr  <- unpacked$cif
  model_name <- unpacked$model_key

  tau <- as.numeric(tau)

  # Snap tau to time_grid once — shared across all causes
  time_grid <- unpacked$time_grid
  snap_idx  <- vapply(tau, function(t) which.min(abs(time_grid - t)),
                      integer(1))

  causes    <- cr@causes
  cause_nms <- paste0("cause_", causes)
  time      <- cr@data[[cr@time_var]]
  status    <- cr@data[[cr@event_var]]
  cens_code <- cr@cens_code

  # Fit marginal model once — shared across all causes and all tau
  margForm <- prodlim::Hist(time, status, cens.code = cens_code) ~ 1
  margFit  <- prodlim::prodlim(margForm, data = cr@data)

  graphs_out          <- stats::setNames(vector("list", length(causes)), cause_nms)
  values_out          <- stats::setNames(vector("list", length(causes)), cause_nms)
  calib_measures_list <- vector("list", length(causes))
  OE_list             <- vector("list", length(causes))

  for (i in seq_along(causes)) {
    k        <- causes[i]
    cause_nm <- cause_nms[i]

    # Slice cause and snap to tau
    pred_full   <- cif_arr[, i, , drop = TRUE]         # [n, Tm_grid]
    pred_at_tau <- pred_full[, snap_idx, drop = FALSE]  # [n, length(tau)]

    res <- .calibration_one_cause(
      pred_at_tau     = pred_at_tau,
      tau             = tau,
      cause           = k,
      cause_idx       = i,
      cause_nm        = cause_nm,
      model_name      = model_name,
      time            = time,
      status          = status,
      cens_code       = cens_code,
      cr              = cr,
      margFit         = margFit,
      loess_smoothing = loess_smoothing,
      bandwidth       = bandwidth,
      graph           = graph
    )

    graphs_out[[cause_nm]]   <- res$graphs
    values_out[[cause_nm]]   <- res$values
    calib_measures_list[[i]] <- res$calib_measures
    OE_list[[i]]             <- res$OE_summary
  }

  list(
    graphs         = if (graph) graphs_out else NULL,
    values         = values_out,
    calib_measures = do.call(rbind, calib_measures_list),
    OE_summary     = do.call(rbind, OE_list)
  )
}
