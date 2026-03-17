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
#'   `"calibration"` (calibration summary statistics),
#'   `"clinical_utility"` (decision curve analysis).
#'   Note: requesting `"IBS"` implicitly also computes `"Brier"`.
#' @param pred_horizons Numeric vector of prediction horizons at which metrics
#'   are evaluated (length `Tm`).  Each value is mapped to the nearest point
#'   on the CIF time grid, so the CIF can be estimated at high resolution
#'   independently of the evaluation horizons.  Required when any of
#'   `"Brier"`, `"IBS"`, `"tdAUC"`, `"cindex_t_year"`, or `"calibration"`
#'   is requested.  Not needed when only `"cindex_rmlt"` is requested.
#' @param args_riskRegression A named list of additional arguments passed to
#'   `riskRegression::Score()`. Relevant for `"Brier"`, `"IBS"`, and `"tdAUC"`.
#'   Defaults: `list(cens.method = "ipcw", cens.model = "km", se.fit = FALSE)`.
#'   Any element supplied here overrides the corresponding default.
#' @param args_pec A named list of additional arguments passed to
#'   `pec::cindex()`. Relevant for `"cindex_t_year"` and `"cindex_rmlt"`.
#'   Defaults: `list(cens.model = "marginal", splitMethod = "noPlan",
#'   verbose = FALSE)`. Any element supplied here overrides the corresponding
#'   default. `eval.times` controls the time horizon used for IPCW weight
#'   truncation; when `NULL`, defaults to the maximum observed event time in
#'   `cr` and a message is issued.
#' @param args_rmlt A named list of additional arguments passed to
#'   [compute_rmlt()]. Relevant for `"cindex_rmlt"`. Defaults:
#'   `list(maxT = NULL)` (resolved to maximum observed event time inside
#'   [compute_rmlt()]). Any element supplied here overrides the corresponding
#'   default.
#' @param args_calibration A named list of arguments passed to the internal
#'   calibration routine. Relevant for `"calibration"`. Defaults:
#'   `list(loess_smoothing = TRUE, bandwidth = NULL, graph = TRUE)`. Any
#'   element supplied here overrides the corresponding default.
#'   `loess_smoothing` controls whether LOESS or nearest-neighbour smoothing
#'   is used; `bandwidth` sets the bandwidth for nearest-neighbour smoothing
#'   (ignored when `loess_smoothing = TRUE`); `graph` controls whether
#'   calibration ggplot objects are produced (set to `FALSE` to skip plot
#'   generation when only `calibration` statistics are needed; `calib_graphs`
#'   will then be absent from the result).
#' @param args_clinical_utility A named list of arguments passed to the
#'   internal DCA routine. Relevant for `"clinical_utility"`. Defaults:
#'   `list(graph = TRUE, xstart = 0.01, xstop = 0.99, xby = 0.01,
#'   ymin = -0.05, harm = 0, intervention = FALSE, interventionper = 100,
#'   smooth = FALSE, loess.span = 0.10)`. Any element supplied here overrides
#'   the corresponding default. `graph` controls whether a ggplot is produced
#'   per cause per horizon (set to `FALSE` to skip plot generation). The
#'   result always contains `dca` — a named list per cause, each a named list
#'   per `pred_horizons` value, each with `net_benefit`,
#'   `interventions_avoided`, and optionally `plot`.
#' @param snap_tol Maximum allowable absolute difference between a
#'   `pred_horizons` value and its nearest point on the CIF time grid.  An
#'   error is raised if any mapped difference exceeds this tolerance (default
#'   `0.1`).  Set to a finite value to guard against accidentally passing
#'   horizons that are far from the estimated time grid, or to `Inf` to
#'   disable the check entirely.
#' @param collapse_as_df Logical; if `TRUE` (default), metrics are collapsed
#'   into data frames:
#'   \describe{
#'     \item{`metrics`}{A data frame with columns `cause`, `pred_horizons`,
#'       and one column per requested horizon-varying metric (`"Brier"`,
#'       `"IBS"`, `"tdAUC"`, `"cindex_t_year"`) and/or calibration statistic
#'       (`ICI`, `E50`, `E90`, `Emax`, `RSB`, `OE`, `OE_lower`, `OE_upper`).
#'       Horizon-varying metrics and calibration columns are joined on
#'       `cause` + `pred_horizons`.}
#'     \item{`scalar_metrics`}{A data frame with columns `cause` and one
#'       column per requested scalar metric (`"cindex_rmlt"`). Only present
#'       if a scalar metric was requested.}
#'     \item{`calib_graphs`}{A named list (keyed by cause) of ggplot objects,
#'       one per `pred_horizons` value. Always a list regardless of
#'       `collapse_as_df`.}
#'   }
#'   When `collapse_as_df = FALSE`, each metric is returned as a named list
#'   with one entry per cause.
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
#'   \item{`"calibration"`}{Calibration summary statistics (ICI, E50, E90,
#'     Emax, RSB) computed via an internal `.compute_calibration()` routine.
#'     These quantify agreement between predicted CIF values and observed
#'     event fractions at each evaluation time.  Calibration plots are
#'     always returned alongside in `calib_graphs`.}
#' }
#'
#' @return When `collapse_as_df = TRUE` (default), a named list with some
#'   subset of the following elements, depending on which metrics were
#'   requested:
#'   \describe{
#'     \item{`metrics`}{Data frame with columns `cause`, `pred_horizons`, and
#'       one column per requested horizon-varying metric and/or calibration
#'       statistic. Present whenever at least one horizon-varying metric or
#'       `"calibration"` was requested.}
#'     \item{`scalar_metrics`}{Data frame with columns `cause` and one column
#'       per requested scalar metric. Present only when `"cindex_rmlt"` is
#'       requested.}
#'     \item{`calib_graphs`}{Named list (keyed by cause) of ggplot objects,
#'       one per `pred_horizons` value. Present only when `"calibration"`
#'       is requested and `args_calibration$graph = TRUE` (the default).}
#'     \item{`dca`}{Named list (keyed by cause), each a named list keyed by
#'       `pred_horizons` value. Each element contains `net_benefit`,
#'       `interventions_avoided`, and optionally `plot`. Present only when
#'       `"clinical_utility"` is requested.}
#'   }
#'   When `collapse_as_df = FALSE`, a named list with one element per
#'   requested metric, each being a named list with one entry per cause.
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
                                                       verbose     = FALSE),
                            args_rmlt           = list(maxT = NULL),
                            args_calibration    = list(loess_smoothing = TRUE,
                                                       bandwidth       = NULL,
                                                       graph           = TRUE),
                            args_clinical_utility = list(graph         = TRUE,
                                                         xstart        = 0.01,
                                                         xstop         = 0.99,
                                                         xby           = 0.01,
                                                         ymin          = -0.05,
                                                         harm          = 0,
                                                         intervention  = FALSE,
                                                         interventionper = 100,
                                                         smooth        = FALSE,
                                                         loess.span    = 0.10),
                            snap_tol       = 0.1,
                            collapse_as_df = TRUE) {
  .check_cr(cr)

  valid_metrics <- c("Brier", "IBS", "tdAUC", "cindex_t_year",
                     "cindex_rmlt", "calibration", "clinical_utility")
  bad_metrics <- setdiff(metrics, valid_metrics)
  if (length(bad_metrics) > 0)
    stop(sprintf(
      "Unrecognised metric(s): %s. Valid options are: %s.",
      paste(bad_metrics, collapse = ", "),
      paste(valid_metrics, collapse = ", ")
    ), call. = FALSE)

  needs_pred_horizons <- any(c("Brier", "IBS", "tdAUC", "cindex_t_year",
                               "calibration", "clinical_utility") %in% metrics)
  if (needs_pred_horizons && is.null(pred_horizons))
    stop(
      '`pred_horizons` must be provided when any of "Brier", "IBS", "tdAUC", ',
      '"cindex_t_year", "calibration", or "clinical_utility" is requested.', call. = FALSE
    )
  if (!is.null(pred_horizons)) {
    if (!is.numeric(pred_horizons) || length(pred_horizons) == 0)
      stop("`pred_horizons` must be a non-empty numeric vector.", call. = FALSE)
    if (anyNA(pred_horizons))
      stop("`pred_horizons` must not contain NA values.", call. = FALSE)
    if (is.unsorted(pred_horizons))
      stop("`pred_horizons` must be sorted in ascending order.", call. = FALSE)
  }

  if (!is.numeric(snap_tol) || length(snap_tol) != 1 || is.na(snap_tol) || snap_tol < 0)
    stop("`snap_tol` must be a single non-negative numeric value.", call. = FALSE)

  valid_calib_args <- c("loess_smoothing", "bandwidth", "graph")
  bad_calib_args   <- setdiff(names(args_calibration), valid_calib_args)
  if (length(bad_calib_args) > 0)
    warning(sprintf(
      "Unrecognised element(s) in `args_calibration` will be ignored: %s. ",
      paste(bad_calib_args, collapse = ", ")
    ), call. = FALSE)
  if (!is.null(args_calibration$loess_smoothing) &&
      (!is.logical(args_calibration$loess_smoothing) ||
       length(args_calibration$loess_smoothing) != 1 ||
       is.na(args_calibration$loess_smoothing)))
    stop("`args_calibration$loess_smoothing` must be a single non-NA logical.",
         call. = FALSE)
  if (!is.null(args_calibration$graph) &&
      (!is.logical(args_calibration$graph) ||
       length(args_calibration$graph) != 1 ||
       is.na(args_calibration$graph)))
    stop("`args_calibration$graph` must be a single non-NA logical.",
         call. = FALSE)
  if (!is.null(args_calibration$bandwidth) &&
      (!is.numeric(args_calibration$bandwidth) ||
       length(args_calibration$bandwidth) != 1 ||
       is.na(args_calibration$bandwidth) ||
       args_calibration$bandwidth <= 0))
    stop("`args_calibration$bandwidth` must be a single positive numeric.",
         call. = FALSE)

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
  rmlt_args  <- modifyList(list(maxT = NULL), args_rmlt)
  calib_args <- modifyList(
    list(loess_smoothing = TRUE, bandwidth = NULL, graph = TRUE),
    args_calibration
  )

  valid_cu_args <- c("graph", "xstart", "xstop", "xby", "ymin", "harm",
                     "intervention", "interventionper", "smooth", "loess.span")
  bad_cu_args <- setdiff(names(args_clinical_utility), valid_cu_args)
  if (length(bad_cu_args) > 0)
    warning(sprintf(
      "Unrecognised element(s) in `args_clinical_utility` will be ignored: %s.",
      paste(bad_cu_args, collapse = ", ")
    ), call. = FALSE)
  cu_args <- modifyList(
    list(graph = TRUE, xstart = 0.01, xstop = 0.99, xby = 0.01,
         ymin = -0.05, harm = 0, intervention = FALSE,
         interventionper = 100, smooth = FALSE, loess.span = 0.10),
    args_clinical_utility
  )

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
  comp_calib     <- "calibration"       %in% metrics
  comp_cu        <- "clinical_utility"  %in% metrics

  if (is.null(pec_args$eval.times) && (comp_cidx_t || comp_cidx_rmlt)) {
    pec_args$eval.times <- max(cr@data[[time_var]][cr@data[[event_var]] != cr@cens_code])
    message(
      "`args_pec$eval.times` not supplied; defaulting to the maximum observed ",
      "event time (", round(pec_args$eval.times, 4), "). This affects C-index ",
      "estimates (`cindex_t_year` and `cindex_rmlt`) only. IPCW weights are ",
      "unstable near the end of follow-up — consider supplying a value via ",
      "`args_pec$eval.times`."
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
  calibration <- named_list(comp_calib)
  calib_graphs   <- named_list(comp_calib)
  dca            <- named_list(comp_cu)

  if (comp_cidx_rmlt) {
    rmlt <- do.call(compute_rmlt, c(
      list(cr  = cr,
           cif = list(cif = cif, time_grid = cif_time_grid)),
      rmlt_args
    ))
  }

  if (!is.null(pred_horizons)) {
    # Warn if any pred_horizons lie outside the CIF time grid range.
    if (any(pred_horizons > max(cif_time_grid)))
      warning(sprintf(
        "The following `pred_horizons` values exceed the maximum CIF time grid point (%g): %s. ",
        max(cif_time_grid),
        paste(pred_horizons[pred_horizons > max(cif_time_grid)], collapse = ", ")
      ), call. = FALSE)
    if (any(pred_horizons < min(cif_time_grid)))
      warning(sprintf(
        "The following `pred_horizons` values are below the minimum CIF time grid point (%g): %s. ",
        min(cif_time_grid),
        paste(pred_horizons[pred_horizons < min(cif_time_grid)], collapse = ", ")
      ), call. = FALSE)

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
  } else {
    snap_idx      <- NULL
    snapped_times <- NULL
  }

  for (k in causes) {
    i  <- which(causes == k)
    nm <- cause_nms[i]

    if (!is.null(snap_idx)) {
      M      <- cif[, i, snap_idx, drop = FALSE]   # [n, length(pred_horizons)]
      dim(M) <- c(nrow(M), length(snap_idx))        # ensure matrix even for n=1
      preds  <- stats::setNames(list(M), nm)
    } else {
      preds <- NULL
    }

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
      brier_sc <- as.data.frame(sc$Brier$score)
      brier_sc <- brier_sc[brier_sc$model == nm & !is.na(brier_sc$times), ]
      if (comp_ibs && is.null(brier_sc$IBS))
        stop(sprintf(
          "Expected an `IBS` column in riskRegression::Score output for cause %s but none found. ",
          "This may indicate a version incompatibility with riskRegression.",
          nm
        ), call. = FALSE)
      if (comp_brier) Brier[[nm]] <- stats::setNames(brier_sc$Brier, pred_horizons)
      if (comp_ibs)   IBS[[nm]]   <- stats::setNames(brier_sc$IBS,   pred_horizons)
      if (comp_auc) {
        auc_sc <- as.data.frame(sc$AUC$score)
        tdAUC[[nm]] <- stats::setNames(
          auc_sc[auc_sc$model == nm, "AUC"],
          pred_horizons
        )
      }
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
      cindex_t_year[[nm]] <- stats::setNames(cidx$AppCindex[[nm]], pred_horizons)
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

    if (comp_cu && !is.null(snap_idx)) {
      ph_nms    <- paste0("pred_horizons_", pred_horizons)
      dca[[nm]] <- stats::setNames(vector("list", length(pred_horizons)), ph_nms)
      for (j in seq_along(pred_horizons)) {
        dca[[nm]][[j]] <- do.call(.clinical_utility, c(
          list(predictions = cif[, i, snap_idx[j]],
               cr          = cr,
               cause       = k,
               cause_nm    = nm,
               timepoint   = snapped_times[j]),
          cu_args
        ))
      }
    }

  }

  if (comp_calib) {
    cal <- .compute_calibration(
      cr              = cr,
      cif_arr         = cif,
      cif_time_grid   = cif_time_grid,
      pred_horizons   = pred_horizons,
      snapped_times   = snapped_times,
      loess_smoothing = calib_args$loess_smoothing,
      bandwidth       = calib_args$bandwidth,
      graph           = calib_args$graph
    )
    for (nm in cause_nms) {
      calibration[[nm]] <- cal$calib_measures[[nm]]
      calib_graphs[[nm]]   <- cal$graphs[[nm]]
    }
  }

  result <- Filter(Negate(is.null), list(
    Brier          = Brier,
    IBS            = IBS,
    tdAUC          = tdAUC,
    cindex_t_year  = cindex_t_year,
    cindex_rmlt    = cindex_rmlt,
    calibration    = calibration,
    dca            = dca
  ))

  # calib_graphs and dca are nested lists — attached at the end
  if (!collapse_as_df) {
    if (comp_calib && calib_args$graph) result$calib_graphs <- calib_graphs
    if (comp_cu) result$dca <- dca
    return(result)
  }

  # dca is a nested list (not collapsible) — extract before collapsing
  dca_result <- result$dca
  result$dca <- NULL

  # All horizon-varying metrics (one value per pred_horizon per cause).
  # Collapsed into a single long data frame: cause, pred_horizons, <metric cols>.
  horizon_metric_nms <- intersect(names(result),
                                  c("Brier", "IBS", "tdAUC", "cindex_t_year"))
  if (length(horizon_metric_nms) > 0) {
    cause_nms_present <- names(result[[horizon_metric_nms[1]]])
    horizon_df <- do.call(rbind, lapply(cause_nms_present, function(cn) {
      row <- data.frame(cause = cn, pred_horizons = pred_horizons)
      for (mn in horizon_metric_nms)
        row[[mn]] <- as.numeric(result[[mn]][[cn]])
      row
    }))
    rownames(horizon_df) <- NULL
    result[horizon_metric_nms] <- NULL
  } else {
    horizon_df <- NULL
  }

  # calibration: named list per cause -> long data frame
  if (!is.null(result$calibration)) {
    lst <- result$calibration
    calib_df <- do.call(rbind, lapply(names(lst), function(nm) {
      cbind(cause = nm, lst[[nm]], row.names = NULL)
    }))
    rownames(calib_df) <- NULL
    result$calibration <- NULL
  } else {
    calib_df <- NULL
  }

  # Merge horizon-varying metrics and calibration into a single data frame,
  # joining on cause + pred_horizons when both are present.
  if (!is.null(horizon_df) && !is.null(calib_df)) {
    result$metrics <- merge(horizon_df, calib_df,
                            by = c("cause", "pred_horizons"), all = TRUE,
                            sort = FALSE)
  } else if (!is.null(horizon_df)) {
    result$metrics <- horizon_df
  } else if (!is.null(calib_df)) {
    result$metrics <- calib_df
  }

  # Scalar metrics (one value per cause): cause, <metric cols>.
  scalar_metric_nms <- intersect(names(result), c("cindex_rmlt"))
  if (length(scalar_metric_nms) > 0) {
    cause_nms_present <- names(result[[scalar_metric_nms[1]]])
    scalar_df <- data.frame(cause = cause_nms_present)
    for (mn in scalar_metric_nms)
      scalar_df[[mn]] <- vapply(result[[mn]], as.numeric, numeric(1))
    rownames(scalar_df) <- NULL
    result[scalar_metric_nms] <- NULL
    result$scalar_metrics <- scalar_df
  }

  # Re-attach calib_graphs and dca unchanged
  if (comp_calib && calib_args$graph) result$calib_graphs <- calib_graphs
  if (comp_cu) result$dca <- dca_result

  result
}
