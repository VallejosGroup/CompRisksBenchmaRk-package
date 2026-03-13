#' Compute performance metrics from CIF predictions
#'
#' Wraps `riskRegression::Score()`, `pec::cindex()`, and custom calibration
#' routines to compute one or more performance metrics for competing risks
#' models. Predictions can be supplied either as a pre-computed CIF array via
#' `cif`, or as a fitted model object via `fit` (in which case the CIF is
#' computed internally using [predict_cif()]).  Exactly one of `cif` or `fit`
#' must be non-`NULL`.
#'
#' @param cr         A [cr_data()] object providing the test-set outcomes, time
#'   and event variable names, and cause codes.
#' @param cif        A list as returned by [predict_cif()], with elements
#'   `$cif` (3-D numeric array `[n, K, Tm]`) and `$time_grid` (numeric vector).
#'   Mutually exclusive with `fit`.
#' @param fit        A fitted model object as returned by [fit_cr_model()].
#'   When supplied, the CIF is computed via [predict_cif()] before scoring.
#'   Mutually exclusive with `cif`.
#' @param cif_time_grid Numeric vector of time points passed to [predict_cif()]
#'   when `fit` is non-`NULL` (default `NULL`).  Must be provided when `fit`
#'   is supplied; must be `NULL` when `cif` is supplied directly (the time grid
#'   is then taken from `cif$time_grid`).
#' @param metrics    Character vector of metrics to compute. Supported values:
#'   `"Brier"` (time-dependent Brier score),
#'   `"IBS"` (integrated Brier score),
#'   `"tdAUC"` (time-dependent AUC),
#'   `"cindex_t_year"` (C-index via \pkg{pec}),
#'   `"cindex_rmlt"` (C-index via restricted mean lifetime),
#'   `"calib_measures"` (calibration summary statistics).
#'   Note: requesting `"IBS"` implicitly also computes `"Brier"`.
#' @param pred_horizons Numeric vector of prediction horizons at which metrics
#'   are evaluated (length `Tm`).  Each value is mapped to the nearest point
#'   on the CIF time grid, so the CIF can be estimated at high resolution
#'   independently of the evaluation horizons.  Required when any of
#'   `"Brier"`, `"IBS"`, `"tdAUC"`, `"cindex_t_year"`, or `"calib_measures"`
#'   is requested.  Not needed when only `"cindex_rmlt"` is requested.
#' @param args_riskRegression A named list of additional arguments passed to
#'   `riskRegression::Score()`. Relevant for `"Brier"`, `"IBS"`, and `"tdAUC"`.
#'   Defaults: `list(cens.method = "ipcw", cens.model = "km", se.fit = FALSE)`.
#'   Any element supplied here overrides the corresponding default.
#' @param args_pec A named list of additional arguments passed to
#'   `pec::cindex()`. Relevant for `"cindex_t_year"` and `"cindex_rmlt"`.
#'   Defaults: `list(cens.model = "marginal", splitMethod = "noPlan",
#'   verbose = FALSE, eval.times = NULL)`. Any element supplied here overrides
#'   the corresponding default.  `eval.times` controls the time horizon used
#'   for IPCW weight truncation in C-index estimation; when `NULL` it defaults
#'   to the maximum observed event time in `cr` and a message is issued, as
#'   IPCW weights are unstable near the end of follow-up.
#' @param args_rmlt A named list of additional arguments passed to
#'   [compute_rmlt()]. Relevant for `"cindex_rmlt"`. Defaults:
#'   `list(maxT = NULL)` (resolved to maximum observed event time inside
#'   [compute_rmlt()]). Any element supplied here overrides the corresponding
#'   default.
#' @param snap_tol Maximum allowable absolute difference between a
#'   `pred_horizons` value and its nearest point on the CIF time grid.  An
#'   error is raised if any mapped difference exceeds this tolerance (default
#'   `0.1`).  Set to a finite value to guard against accidentally passing
#'   horizons that are far from the estimated time grid, or to `Inf` to
#'   disable the check entirely.
#' @param collapse_as_df Logical; if `TRUE` (default), each metric element in
#'   the output list is collapsed into a data frame with one row per cause.
#'   For time-varying metrics (`"Brier"`, `"tdAUC"`, `"cindex_t_year"`),
#'   columns correspond to `pred_horizons`.  For scalar metrics (`"IBS"`,
#'   `"cindex_rmlt"`), a single `value` column is used.  For
#'   `"calib_measures"`, the result is a data frame with columns `cause`,
#'   `pred_horizons`, and the calibration statistics (`ICI`, `E50`, `E90`,
#'   `Emax`, `RSB`), with one row per cause per prediction horizon.
#'   `calib_graphs` is always returned as a named list regardless of this
#'   setting.
#'
#' @details
#' The following metrics are supported, computed separately for each competing
#' cause:
#'
#' \describe{
#'   \item{`"Brier"`}{The time-dependent Brier score, an IPCW-weighted mean
#'     squared error between predicted CIF values and event indicators at each
#'     evaluation time.  Lower values indicate better calibration and
#'     discrimination.  Computed via `riskRegression::Score()`.}
#'
#'   \item{`"IBS"`}{The integrated Brier score, summarising the time-dependent
#'     Brier score into a single value by integration over `pred_horizons`.
#'     Requesting `"IBS"` implicitly triggers computation of `"Brier"` as well.
#'     Computed via `riskRegression::Score()`.}
#'
#'   \item{`"tdAUC"`}{The time-dependent cause-specific AUC, measuring the
#'     probability that a randomly chosen case has a higher predicted CIF than
#'     a randomly chosen control at each evaluation time.  Values above 0.5
#'     indicate discrimination better than chance.  Computed via
#'     `riskRegression::Score()`.}
#'
#'   \item{`"cindex_t_year"`}{The concordance index (C-index) computed via
#'     `pec::cindex()` using a marginal censoring model.  Summarises
#'     discrimination across the entire follow-up rather than at a specific
#'     time point.}
#'
#'   \item{`"cindex_rmlt"`}{A C-index based on the restricted mean lifetime
#'     (RMLT), computed by first integrating each subject's CIF up to
#'     `args_rmlt$maxT` via trapezoidal integration, then passing the resulting
#'     summary scores to `pec::cindex()`.  Useful as a scalar discrimination
#'     summary that does not depend on a specific time horizon.}
#'
#'   \item{`"calib_measures"`}{Calibration summary statistics (ICI, E50, E90,
#'     Emax, RSB) computed via an internal `.compute_calibration()` routine.
#'     These quantify agreement between predicted CIF values and observed
#'     event fractions at each evaluation time.  Calibration plots are
#'     always returned alongside in `calib_graphs`.}
#' }
#'
#' @return A named list with one element per requested metric.  When
#'   `collapse_as_df = TRUE` (default), each element (except `calib_graphs`)
#'   is a data frame with one row per cause.  When `collapse_as_df = FALSE`,
#'   each element is itself a named list with one entry per cause.
#'   `calib_graphs` is always a named list keyed by cause (e.g. `"cause_1"`),
#'   each containing a list of ggplot objects — one per element of
#'   `pred_horizons`.
#' @export
compute_metrics <- function(cr,
                            cif           = NULL,
                            fit           = NULL,
                            cif_time_grid = NULL,
                            metrics       = c("Brier", "IBS", "tdAUC", "cindex_rmlt"),
                            pred_horizons = NULL,
                            args_riskRegression = list(cens.method = "ipcw",
                                                       cens.model  = "km",
                                                       se.fit      = FALSE),
                            args_pec            = list(cens.model  = "marginal",
                                                       splitMethod = "noPlan",
                                                       verbose     = FALSE,
                                                       eval.times  = NULL),
                            args_rmlt           = list(maxT = NULL),
                            snap_tol       = 0.1,
                            collapse_as_df = TRUE) {
  .check_cr(cr)

  valid_metrics <- c("Brier", "IBS", "tdAUC", "cindex_t_year",
                     "cindex_rmlt", "calib_measures")
  bad_metrics <- setdiff(metrics, valid_metrics)
  if (length(bad_metrics) > 0)
    stop(sprintf(
      "Unrecognised metric(s): %s. Valid options are: %s.",
      paste(bad_metrics, collapse = ", "),
      paste(valid_metrics, collapse = ", ")
    ), call. = FALSE)

  needs_pred_horizons <- any(c("Brier", "IBS", "tdAUC", "cindex_t_year",
                             "calib_measures") %in% metrics)
  if (needs_pred_horizons && is.null(pred_horizons))
    stop(
      '`pred_horizons` must be provided when any of "Brier", "IBS", "tdAUC", ',
      '"cindex_t_year", or "calib_measures" is requested.', call. = FALSE
    )
  if (!is.null(pred_horizons)) {
    if (!is.numeric(pred_horizons) || length(pred_horizons) == 0)
      stop("`pred_horizons` must be a non-empty numeric vector.", call. = FALSE)
    if (anyNA(pred_horizons))
      stop("`pred_horizons` must not contain NA values.", call. = FALSE)
    if (is.unsorted(pred_horizons))
      stop("`pred_horizons` must be sorted in ascending order.", call. = FALSE)
  }

  cif               <- .resolve_cif(cif, fit, cif_time_grid, cr)
  unpacked          <- .validate_and_unpack_cif(cif, cr)
  cif_time_grid     <- unpacked$time_grid
  cif               <- unpacked$cif

  # Merge user-supplied args over defaults
  rr_args   <- modifyList(
    list(cens.method = "ipcw", cens.model = "km", se.fit = FALSE),
    args_riskRegression
  )
  pec_args  <- modifyList(
    list(cens.model = "marginal", splitMethod = "noPlan", verbose = FALSE),
    args_pec
  )
  rmlt_args <- modifyList(list(maxT = NULL), args_rmlt)

  time_var  <- cr@time_var
  event_var <- cr@event_var
  causes    <- cr@causes
  cause_nms <- paste0("cause_", causes)

  # Build a Hist() formula for pec::cindex(); environment must be set to
  # prodlim so that Hist() resolves correctly inside riskRegression/pec calls.
  f <- stats::as.formula(sprintf("Hist(%s,%s) ~ 1", time_var, event_var))
  environment(f) <- asNamespace("prodlim")

  comp_brier     <- "Brier"          %in% metrics
  comp_ibs       <- "IBS"            %in% metrics
  comp_auc       <- "tdAUC"          %in% metrics
  comp_cidx_t    <- "cindex_t_year"  %in% metrics
  if (comp_cidx_t)
    message(
      "`cindex_t_year` is not recommended: the C-index is not a proper metric ",
      "for t-year predictions. Consider `cindex_rmlt` instead."
    )
  comp_cidx_rmlt <- "cindex_rmlt"    %in% metrics
  comp_calib     <- "calib_measures" %in% metrics

  if (is.null(pec_args$eval.times) && (comp_cidx_t || comp_cidx_rmlt)) {
    pec_args$eval.times <- max(cr@data[[time_var]][cr@data[[event_var]] != cr@cens_code])
    message(
      "`args_pec$eval.times` not supplied; defaulting to the maximum observed ",
      "event time (", round(pec_args$eval.times, 4), "). This affects C-index ",
      "estimates (`cindex_t_year` and `cindex_rmlt`) only. IPCW weights are ",
      "unstable near the end of follow-up — consider supplying a lower value."
    )
  }

  # IBS requires Brier from riskRegression; always request "brier" when either
  # "Brier" or "IBS" is asked for. Translate user-facing names to riskRegression
  # expected strings.
  need_rr_brier <- comp_brier || comp_ibs
  rr_metrics <- c(
    if (need_rr_brier) "brier",
    if (comp_auc)      "auc"
  )
  rr_summary <- if (comp_ibs) "ibs" else character(0)

  named_list <- function(cond)
    if (cond) stats::setNames(vector("list", length(causes)), cause_nms) else NULL

  Brier          <- named_list(comp_brier)
  IBS            <- named_list(comp_ibs)
  tdAUC          <- named_list(comp_auc)
  cindex_t_year  <- named_list(comp_cidx_t)
  cindex_rmlt    <- named_list(comp_cidx_rmlt)
  calib_measures <- named_list(comp_calib)
  calib_graphs   <- named_list(comp_calib)

  if (comp_cidx_rmlt) {
    rmlt <- do.call(compute_rmlt, c(
      list(cr  = cr,
           cif = list(cif = cif, time_grid = cif_time_grid)),
      rmlt_args
    ))
  }

  # Map each pred_horizons value to the nearest point on the CIF time grid.
  # riskRegression::Score and pec::cindex expect predictions whose columns
  # correspond exactly to the evaluation times supplied; by subsetting the CIF
  # array to the snapped indices we decouple the CIF time resolution from the
  # evaluation times requested by the user.
  snap_idx      <- vapply(pred_horizons,
                          function(t) which.min(abs(cif_time_grid - t)),
                          integer(1))
  snapped_times <- cif_time_grid[snap_idx]
  snap_diffs    <- abs(pred_horizons - snapped_times)
  if (any(snap_diffs > snap_tol))
    stop(sprintf(
      paste0("The following `pred_horizons` values are more than `snap_tol` = %g ",
             "away from the nearest CIF time grid point:\n%s\n",
             "Adjust `pred_horizons`, increase `snap_tol`, or use a finer CIF time grid."),
      snap_tol,
      paste(sprintf("  %g (nearest grid point: %g, diff: %g)",
                    pred_horizons[snap_diffs > snap_tol],
                    snapped_times[snap_diffs > snap_tol],
                    snap_diffs[snap_diffs > snap_tol]),
            collapse = "\n")
    ), call. = FALSE)

  for (k in causes) {
    i      <- which(causes == k)
    nm     <- cause_nms[i]
    M      <- cif[, i, snap_idx, drop = FALSE]   # [n, length(pred_horizons)]
    dim(M) <- c(nrow(M), length(snap_idx))        # ensure matrix even for n=1
    preds  <- stats::setNames(list(M), nm)

    if (length(rr_metrics) > 0) {
      sc <- do.call(riskRegression::Score, c(
        list(object  = preds,
             formula = f,
             data    = cr@data,
             metrics = rr_metrics,
             summary = rr_summary,
             times   = snapped_times,
             cause   = k,
             null.model = FALSE),
        rr_args
      ))
      if (comp_brier) Brier[[nm]] <- sc$Brier$score[sc$Brier$score$model == nm, "Brier"]
      if (comp_ibs)   IBS[[nm]]   <- sc$Brier$score[sc$Brier$score$model == nm, "IBS"]
      if (comp_auc)   tdAUC[[nm]] <- sc$AUC$score[sc$AUC$score$model == nm,     "AUC"]
    }

    if (comp_cidx_t) {
      pec_args_cidx_t <- modifyList(
        pec_args,
        list(eval.times = rep(pec_args$eval.times, length(snapped_times)),
             pred.times = snapped_times)
      )
      cidx <- do.call(pec::cindex, c(
        list(object  = preds,
             formula = f,
             data    = cr@data,
             cause   = k),
        pec_args_cidx_t
      ))
      cindex_t_year[[nm]] <- cidx$AppCindex[[nm]]
    }
    if (comp_cidx_rmlt) {
      preds_rmlt           <- as.matrix(rmlt[, i])
      colnames(preds_rmlt) <- nm
      preds_rmlt           <- stats::setNames(list(preds_rmlt), nm)
      cidx_rmlt <- do.call(pec::cindex, c(
        list(object     = preds_rmlt,
             formula    = f,
             data       = cr@data,
             cause      = k),
        pec_args
      ))
      cindex_rmlt[[nm]] <- cidx_rmlt$AppCindex[[nm]]
    }

  }

  if (comp_calib) {
    cal <- .compute_calibration(
      cr              = cr,
      cif_arr         = cif,
      cif_time_grid   = cif_time_grid,
      pred_horizons   = snapped_times,
      loess_smoothing = TRUE,
      graph           = TRUE
    )
    for (nm in cause_nms) {
      calib_measures[[nm]] <- cal$calib_measures[[nm]]
      calib_graphs[[nm]]   <- cal$graphs[[nm]]
    }
  }

  result <- Filter(Negate(is.null), list(
    Brier          = Brier,
    IBS            = IBS,
    tdAUC          = tdAUC,
    cindex_t_year  = cindex_t_year,
    cindex_rmlt    = cindex_rmlt,
    calib_measures = calib_measures
  ))

  # calib_graphs is always a named list (not collapsible) — kept separate
  if (comp_calib) result$calib_graphs <- calib_graphs

  if (!collapse_as_df) return(result)

  t_nms <- if (!is.null(pred_horizons)) paste0("pred_horizons_", seq_along(pred_horizons)) else character(0)

  # Collapse all metrics except calib_graphs (plots cannot be data frames)
  collapse_nms <- setdiff(names(result), "calib_graphs")
  collapsed <- lapply(collapse_nms, function(metric_nm) {
    lst <- result[[metric_nm]]
    if (metric_nm == "cindex_rmlt") {
      vals <- vapply(lst, as.numeric, numeric(1))
      df   <- data.frame(value = vals, row.names = names(lst))
    } else if (metric_nm == "calib_measures") {
      df <- do.call(rbind, lapply(names(lst), function(nm) {
        cbind(cause = nm, lst[[nm]])
      }))
      rownames(df) <- NULL
    } else {
      mat <- do.call(rbind, lapply(lst, as.numeric))
      rownames(mat) <- names(lst)
      df <- as.data.frame(mat)
      colnames(df) <- t_nms
    }
    df
  }) |> stats::setNames(collapse_nms)

  # Re-attach calib_graphs unchanged
  if (comp_calib) collapsed$calib_graphs <- calib_graphs
  collapsed
}
