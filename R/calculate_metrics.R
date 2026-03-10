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
#' @param eval_times Numeric vector of evaluation times (length `Tm`).
#' @param cif        3-D numeric array with dimensions `[n, K, Tm]`, where `n`
#'   is the number of subjects, `K` the number of competing causes, and `Tm`
#'   the number of evaluation times.  Mutually exclusive with `fit`.
#' @param fit        A fitted model object as returned by [fit_cr_model()] or a
#'   registered model's `fit` function.  When supplied, the CIF is computed via
#'   [predict_cif()] before scoring.  Mutually exclusive with `cif`.
#' @param cif_time_grid Numeric vector of time points passed to [predict_cif()]
#'   when `fit` is non-`NULL` (default `NULL`).  Must be provided when `fit`
#'   is supplied; ignored when `cif` is supplied directly.
#' @param metrics    Character vector of metrics to compute. Supported values:
#'   `"brier"` (Brier score), `"auc"` (time-dependent AUC), `"cidx_pec"`
#'   (C-index via \pkg{pec}), `"cidx_survM"` (C-index via SurvM),
#'   `"cidx_rmlt"` (C-index via restricted mean lifetime),
#'   `"calib_measures"` (calibration summary statistics).
#' @param summary    Character vector of Brier-score summaries (relevant only
#'   when `"brier"` is in `metrics`): `"ibs"` for the integrated Brier score
#'   and/or `"risks"` for per-risk summaries.
#' @param cens.method Censoring correction method passed to
#'   `riskRegression::Score()` (default `"ipcw"`); relevant for `"brier"`
#'   and `"auc"`.
#' @param cens.model  Censoring model used for IPCW weighting (default
#'   `"km"`); relevant for `"brier"` and `"auc"`.
#' @param cens.code   Integer code denoting censored observations (default
#'   `0`); relevant for `"calib_measures"`.
#' @param se.fit      Logical; compute standard errors? (default `FALSE`);
#'   relevant for `"brier"` and `"auc"`.
#' @param rmlt_tau    Upper integration limit for the restricted mean lifetime
#'   (default `NULL`, which uses `max(eval_times)`); relevant only for
#'   `"cidx_rmlt"`.
#'
#' @details
#' The following metrics are supported, computed separately for each competing
#' cause:
#'
#' \describe{
#'   \item{`"brier"`}{The time-dependent Brier score, an IPCW-weighted mean
#'     squared error between predicted CIF values and event indicators at each
#'     evaluation time.  Lower values indicate better calibration and
#'     discrimination.  Computed via `riskRegression::Score()`.}
#'
#'   \item{`"auc"`}{The time-dependent cause-specific AUC, measuring the
#'     probability that a randomly chosen case has a higher predicted CIF than
#'     a randomly chosen control at each evaluation time.  Values above 0.5
#'     indicate discrimination better than chance.  Computed via
#'     `riskRegression::Score()`.}
#'
#'   \item{`"cidx_pec"`}{The concordance index (C-index) computed via
#'     `pec::cindex()` using a marginal censoring model.  Summarises
#'     discrimination across the entire follow-up rather than at a specific
#'     time point.}
#'
#'   \item{`"cidx_survM"`}{A cause-specific C-index computed via an internal
#'     `CindexCR()` routine.  Evaluated at each time in `eval_times` and
#'     returned together with the corresponding evaluation time points.}
#'
#'   \item{`"cidx_rmlt"`}{A C-index based on the restricted mean lifetime
#'     (RMLT), computed by first integrating each subject's CIF up to
#'     `rmlt_tau` via trapezoidal integration, then passing the resulting
#'     summary scores to `pec::cindex()`.  Useful as a scalar discrimination
#'     summary that does not depend on a specific time horizon.}
#'
#'   \item{`"calib_measures"`}{Calibration summary statistics (ICI, E50, E90,
#'     Emax, RSB) computed via an internal `CalibrationPlot()` routine.
#'     These quantify agreement between predicted CIF values and observed
#'     event fractions at each evaluation time.}
#' }
#'
#' When `"brier"` is requested, the `summary` argument controls whether the
#' integrated Brier score (`"ibs"`) and/or per-risk summaries (`"risks"`) are
#' also returned.
#'
#' @return A named list; elements present depend on the requested metrics:
#'   `Brier`, `IBS`, `tdAUC`, `cindex_t_year`, `cindex_survM`, `cindex_rmlt`,
#'   `calib_measures`.  Each element is itself a named list with one entry per
#'   cause (e.g. `$Brier$cause_1`, `$Brier$cause_2`).
#' @export
calculate_metrics <- function(cr, eval_times,
                              cif          = NULL,
                              fit          = NULL,
                              cif_time_grid = NULL,
                              metrics      = c("brier", "auc"),
                              summary      = c("ibs", "risks"),
                              cens.method  = "ipcw",
                              cens.model   = "km",
                              cens.code    = 0,
                              se.fit       = FALSE,
                              rmlt_tau     = NULL) {
  if (!methods::is(cr, "cr_data"))
    stop("`cr` must be a cr_data object.", call. = FALSE)

  if (is.null(cif) && is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are non-NULL.",
         call. = FALSE)

  if (!is.null(fit)) {
    if (is.null(cif_time_grid))
      stop("`cif_time_grid` must be provided when `fit` is non-NULL.", call. = FALSE)
    cif <- predict_cif(fit, newdata = cr, time_grid = cif_time_grid)
  }

  time_var  <- cr@time_var
  event_var <- cr@event_var
  causes    <- cr@causes
  idx_nms   <- paste0("cause_", causes)

  d  <- dim(cif)
  n  <- d[1]; Tm <- d[3]

  f <- stats::as.formula(
    sprintf("Hist(%s,%s) ~ 1", time_var, event_var),
    env = getNamespace("prodlim")
  )

  comp_brier <- "brier"          %in% metrics
  comp_auc   <- "auc"            %in% metrics
  comp_pec   <- "cidx_pec"       %in% metrics
  comp_survM <- "cidx_survM"     %in% metrics
  comp_rmlt  <- "cidx_rmlt"      %in% metrics
  comp_calib <- "calib_measures" %in% metrics
  comp_ibs   <- comp_brier && ("ibs" %in% summary)

  rr_metrics <- intersect(metrics, c("brier", "auc"))
  rr_summary <- intersect(summary, c("ibs", "risks"))

  named_list <- function(cond)
    if (cond) stats::setNames(vector("list", length(causes)), idx_nms) else NULL

  Brier          <- named_list(comp_brier)
  IBS            <- named_list(comp_ibs)
  tdAUC          <- named_list(comp_auc)
  cindex_t_year  <- named_list(comp_pec)
  cindex_survM   <- named_list(comp_survM)
  cindex_rmlt    <- named_list(comp_rmlt)
  calib_measures <- named_list(comp_calib)

  if (comp_rmlt) {
    tau  <- if (!is.null(rmlt_tau)) rmlt_tau else max(eval_times)
    rmlt <- compute_rmlt_from_cif(cif, eval_times, tau = tau)
  }

  for (k in causes) {
    i      <- which(causes == k)
    nm     <- idx_nms[i]
    M      <- cif[, i, , drop = FALSE]
    dim(M) <- c(n, Tm)
    k_name <- paste0("cif_cause_", k)
    preds  <- stats::setNames(list(M), k_name)

    if (length(rr_metrics) > 0) {
      sc   <- riskRegression::Score(
        object      = preds,
        formula     = f,
        data        = cr@data,
        metrics     = rr_metrics,
        summary     = rr_summary,
        times       = eval_times,
        cause       = k,
        cens.method = cens.method,
        cens.model  = cens.model,
        se.fit      = se.fit
      )
      rows <- sc$Brier$score$model == k_name
      if (comp_brier) Brier[[nm]]  <- sc$Brier$score[rows, "Brier"]
      if (comp_ibs)   IBS[[nm]] <- sc$Brier$score[rows, "IBS"]
      if (comp_auc)   tdAUC[[nm]] <- sc$AUC$score[
        sc$AUC$score$model == k_name, "AUC"]
    }

    if (comp_pec) {
      cidx <- pec::cindex(
        object      = preds,
        formula     = f,
        data        = cr@data,
        eval.times  = eval_times,
        pred.times  = eval_times,
        cause       = k,
        cens.model  = "marginal",
        splitMethod = "noPlan",
        verbose     = FALSE
      )
      cindex_t_year[[nm]] <- cidx$AppCindex[[k_name]]
    }

    if (comp_survM) {
      res <- lapply(seq_len(ncol(M)), function(j)
        CindexCR(time      = cr@data[[time_var]],
                 status    = cr@data[[event_var]],
                 predicted = 1 - M[, j],
                 Cause_int = k)
      )
      cindex_survM[[nm]] <- list(
        sapply(res, `[[`, "cindex"),
        sapply(res, `[[`, "time_ev")
      )
    }

    if (comp_rmlt) {
      k_name_rmlt <- paste0("cause_", k)
      pred_mat    <- as.matrix(rmlt[, i])
      colnames(pred_mat) <- paste0("times_", max(eval_times))
      preds_rmlt  <- stats::setNames(list(pred_mat), k_name_rmlt)
      cidx_rmlt   <- pec::cindex(
        object      = preds_rmlt,
        formula     = f,
        data        = cr@data,
        eval.times  = max(eval_times),
        cens.model  = "marginal",
        splitMethod = "noPlan",
        verbose     = FALSE
      )
      cindex_rmlt[[nm]] <- cidx_rmlt$AppCindex[[k_name_rmlt]]
    }

    if (comp_calib) {
      cm <- CalibrationPlot(
        predictions      = M,
        data             = cr@data,
        time             = cr@data[[time_var]],
        status           = cr@data[[event_var]],
        tau              = eval_times,
        cause            = k,
        cens.code        = cens.code,
        predictions.type = "CIF",
        loess_smoothing  = TRUE,
        graph            = FALSE
      )$calib.measures
      calib_measures[[nm]] <- do.call(rbind, Filter(Negate(is.null), cm))
    }
  }

  Filter(Negate(is.null), list(
    Brier          = Brier,
    IBS            = IBS,
    tdAUC          = tdAUC,
    cindex_t_year  = cindex_t_year,
    cindex_survM   = cindex_survM,
    cindex_rmlt    = cindex_rmlt,
    calib_measures = calib_measures
  ))
}
