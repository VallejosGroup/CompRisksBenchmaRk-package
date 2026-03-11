#' @title Performance Metrics for Competing Risks Models
#' @description Functions for computing discrimination, calibration, and
#'   clinical utility metrics for competing risks survival models.
#' @name metrics
NULL




#' Weighted Brier Score for competing risks
#'
#' Computes the inverse probability of censoring weighted (IPCW) Brier
#' score at one or more evaluation times.
#'
#' @param predictions Numeric matrix `[n, length(tau)]` of predicted
#'   cumulative incidences (or a numeric vector for a single time point).
#' @param tau  Numeric scalar or vector of evaluation times.
#' @param time Numeric vector of observed times.
#' @param status Integer vector of event codes.
#' @param cause Integer code of the event of interest.
#' @param cens.code Integer code for censoring (default `0`).
#' @param cmprsk Logical; if `TRUE`, use competing risks Brier score
#'   weights (default `FALSE`).
#'
#' @return A list with elements `weighted.brier.score`, `tau`, `n`, and
#'   `n.risk`.
#' @export
WeightedBrierScore <- function(predictions,
                                tau,
                                time,
                                status,
                                cause,
                                cens.code = 0L,
                                cmprsk    = FALSE) {
  if (is.vector(predictions))
    predictions <- matrix(as.numeric(predictions), ncol = 1L)
  else
    predictions <- as.matrix(predictions)
  tau <- as.numeric(tau)

  if (ncol(predictions) == 1L && length(tau) > 1L)
    predictions <- matrix(predictions[, 1L],
                           nrow = nrow(predictions), ncol = length(tau))
  if (ncol(predictions) != length(tau))
    stop("ncol(predictions) must equal length(tau).")

  G <- censor.prob.KM(time = time, status = status, cens.code = cens.code)
  n <- nrow(predictions)
  BS_we_all  <- numeric(length(tau))
  n_risk_all <- integer(length(tau))

  for (jj in seq_along(tau)) {
    tau_j     <- tau[jj]
    p_j       <- predictions[, jj]
    residuals <- rep(0, n)

    for (i in seq_len(n)) {
      indx1 <- which(G[, 1L] >= time[i])
      indx2 <- which(G[, 1L] >= tau_j)
      G1    <- if (length(indx1) > 0L) G[indx1[1L], 2L] else 1
      G2    <- if (length(indx2) > 0L) G[indx2[1L], 2L] else 1

      if (cmprsk) {
        if      (status[i] == cause && time[i] <= tau_j)
          residuals[i] <- (1 - p_j[i])^2 / G1
        else if (status[i] != cause && status[i] != cens.code && time[i] <= tau_j)
          residuals[i] <- p_j[i]^2 / G1
        else if (status[i] == cens.code && time[i] <= tau_j)
          residuals[i] <- 0
        else
          residuals[i] <- p_j[i]^2 / G2
      } else {
        if      (time[i] <= tau_j && status[i] == cause)
          residuals[i] <- (1 - p_j[i])^2 / G1
        else if (time[i] <= tau_j && status[i] == cens.code)
          residuals[i] <- 0
        else
          residuals[i] <- p_j[i]^2 / G2
      }
    }
    BS_we_all[jj]  <- mean(residuals)
    n_risk_all[jj] <- sum(time > tau_j)
  }
  list(weighted.brier.score = BS_we_all, tau = tau, n = n,
       n.risk = n_risk_all)
}


#' Integrated Brier Score for competing risks
#'
#' Wraps [WeightedBrierScore()] and integrates across time using the
#' trapezoidal rule.
#'
#' @param prediction_matrix Numeric matrix `[n, Tm]` of predictions.
#' @param taus   Numeric vector of evaluation times (length `Tm`).
#' @param time   Numeric vector of observed times.
#' @param status Integer vector of event codes.
#' @param cause  Integer code of the event of interest.
#' @param cens.code Integer code for censoring.
#' @param cmprsk Logical; use competing risks Brier score?
#'
#' @return A list with `integrated.brier.score`, `average.brier.score`,
#'   `weighted.brier.score`, and `taus`.
#' @export
IntegratedBrierScore <- function(prediction_matrix,
                                  taus,
                                  time,
                                  status,
                                  cause,
                                  cens.code,
                                  cmprsk) {
  bs <- WeightedBrierScore(
    predictions = prediction_matrix,
    tau         = taus,
    time        = time,
    status      = status,
    cause       = cause,
    cens.code   = cens.code,
    cmprsk      = cmprsk
  )$weighted.brier.score

  t_max <- max(taus)
  t_min <- min(taus)
  ibs   <- .trapezoidal_integration(taus, bs) / (t_max - t_min)

  list(integrated.brier.score = ibs,
       average.brier.score    = mean(bs),
       weighted.brier.score   = bs,
       taus                   = taus)
}


#' @noRd
censor.prob.KM <- function(time, status, cens.code) {
  tmp              <- data.frame(time = time)
  tmp$censor.status <- ifelse(status == cens.code, 0L, 1L)
  fit <- prodlim::prodlim(
    formula = survival::Surv(time, censor.status) ~ 1,
    data    = tmp,
    reverse = TRUE
  )
  prob <- stats::predict(fit,
                         times       = sort(unique(time)),
                         level.chaos = 1,
                         mode        = "matrix",
                         type        = "surv")
  out <- cbind(sort(unique(time)), prob)
  stats::na.omit(out)
}


#' C-index for competing risks (Wolbers/SurvMetrics method)
#'
#' @param time       Numeric vector of observed times.
#' @param status     Integer vector of event codes (0 = censored).
#' @param predicted  Numeric vector of predicted CIF values (or
#'   `1 - CIF` depending on convention; see `Cause_int`).
#' @param Cause_int  Integer indicating the cause of interest (default 1).
#'
#' @return A list with `cindex` (numeric) and `time_ev` (the evaluation
#'   time used internally, i.e. `max(time) + 1`).
#' @export
CindexCR <- function(time, status, predicted, Cause_int = 1L) {
  if (any(is.na(time)))      stop("time cannot contain NA.")
  if (any(is.na(status)))    stop("status cannot contain NA.")
  if (any(!(status %in% c(0L, 1L, 2L))))
    stop("status must be 0, 1, or 2.")
  if (any(is.na(predicted))) stop("predicted cannot contain NA.")
  if (!(Cause_int %in% status)) stop("Invalid Cause_int.")
  if (min(time) <= 0) stop("Survival times must be positive.")

  Censoring  <- ifelse(status == 0L, 0L, 1L)
  Cause      <- ifelse(status == 2L, 2L, 1L)
  Prediction <- -log(predicted)
  Time       <- max(time) + 1

  n          <- length(Prediction)
  A  <- B  <- Q  <- N_t  <- matrix(0L, n, n)

  for (i in seq_len(n)) {
    A[i, which(time[i] < time)] <- 1L
    B[i, intersect(
      intersect(which(time[i] >= time), which(Cause != Cause_int)),
      which(Censoring == 1L)
    )] <- 1L
    Q[i, which(Prediction[i] > Prediction)] <- 1L
    if (time[i] <= Time && Cause[i] == Cause_int && Censoring[i] == 1L)
      N_t[i, ] <- 1L
  }

  Num <- sum((A + B) * Q * N_t)
  Den <- sum((A + B) * N_t)
  list(cindex = Num / Den, time_ev = Time)
}


#' Calibration plot for competing risks
#'
#' Produces a pseudo-value calibration plot using LOESS or nearest-neighbour
#' smoothing.  Returns the plot objects, smoothed calibration frames, and
#' numerical calibration measures.
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
#' @param cause         Integer cause code of interest (must be in
#'   `cr@causes`).
#' @param loess_smoothing If `TRUE` (default), use LOESS smoothing;
#'   otherwise use `prodlim::meanNeighbors()`.
#' @param bandwidth     Optional bandwidth for nearest-neighbour smoothing.
#' @param graph         Logical; produce ggplot objects? (default `TRUE`).
#'
#' @return A list with elements `graphs`, `values`, `calib.measures`, and
#'   `OE_summary`.
#' @export
CalibrationPlot <- function(cr,
                             cif           = NULL,
                             fit           = NULL,
                             cif_time_grid = NULL,
                             tau,
                             cause,
                             loess_smoothing = TRUE,
                             bandwidth       = NULL,
                             graph           = TRUE) {
  .check_cr(cr)

  if (missing(tau) || !is.numeric(tau) || length(tau) == 0)
    stop("`tau` must be a non-empty numeric vector.", call. = FALSE)
  if (anyNA(tau) || is.unsorted(tau))
    stop("`tau` must be sorted and contain no NA values.", call. = FALSE)

  if (missing(cause) || !is.numeric(cause) || length(cause) != 1 || is.na(cause))
    stop("`cause` must be a single numeric cause code.", call. = FALSE)
  cause <- as.integer(cause)
  if (!cause %in% cr@causes)
    stop(sprintf("`cause` %d is not in cr@causes (%s).",
                 cause, paste(cr@causes, collapse = ", ")), call. = FALSE)

  cif      <- .resolve_cif(cif, fit, cif_time_grid, cr)
  unpacked <- .validate_and_unpack_cif(cif, cr)
  cif_arr  <- unpacked$cif
  model_name <- unpacked$model_key

  # Slice the cause dimension
  cause_idx <- which(cr@causes == cause)
  pred_mat  <- cif_arr[, cause_idx, , drop = TRUE]  # [n, Tm]
  if (is.null(dim(pred_mat)))
    pred_mat <- matrix(pred_mat, nrow = nrow(cif_arr))

  tau <- as.numeric(tau)
  if (ncol(pred_mat) != length(unpacked$time_grid))
    stop("CIF time grid length does not match predictions.", call. = FALSE)

  # Interpolate predictions to requested tau from cif time_grid
  time_grid <- unpacked$time_grid
  pred_at_tau <- matrix(NA_real_, nrow = nrow(pred_mat), ncol = length(tau))
  for (j in seq_along(tau)) {
    idx <- which.min(abs(time_grid - tau[j]))
    pred_at_tau[, j] <- pred_mat[, idx]
  }

  time_var  <- cr@time_var
  event_var <- cr@event_var
  time      <- cr@data[[time_var]]
  status    <- cr@data[[event_var]]
  cens_code <- cr@cens_code

  # O/E at the last (or only) horizon
  horizon     <- tau[length(tau)]
  obj         <- summary(
    survival::survfit(
      survival::Surv(time, factor(status, levels = c(cens_code, cr@causes))) ~ 1
    ),
    times = horizon
  )
  aj          <- list(obs = obj$pstate[, cause_idx + 1L],
                      se  = obj$std.err[, cause_idx + 1L])
  horizon_idx <- which.min(abs(tau - horizon))
  OE          <- aj$obs / mean(pred_at_tau[, horizon_idx])
  OE_summary  <- c(
    OE    = OE,
    lower = exp(log(OE) - stats::qnorm(0.975) * aj$se / aj$obs),
    upper = exp(log(OE) + stats::qnorm(0.975) * aj$se / aj$obs)
  )
  e_txt <- with(as.list(OE_summary),
                sprintf(", OE = %.3f (CI: %.3f\u2013%.3f)", OE, lower, upper))

  plotFrames <- vector("list", length(tau))
  measures   <- vector("list", length(tau))
  graphs     <- vector("list", length(tau))

  # Fit marginal model once — loop-invariant
  margForm <- prodlim::Hist(time, status, cens.code = cens_code) ~ 1
  margFit  <- prodlim::prodlim(margForm, data = cr@data)

  for (i in seq_along(tau)) {
    prediction <- pred_at_tau[, i]
    eval_time  <- tau[i]

    if (length(unique(stats::na.omit(prediction))) <= 1) {
      measures[[i]] <- data.frame(tau = eval_time, ICI = NA, E50 = NA,
                                   E90 = NA, Emax = NA, RSB = NA)
      next
    }

    pseudo <- prodlim::jackknife(margFit, cause = cause, times = eval_time)

    keep <- !is.na(prediction) & !is.na(pseudo)
    if (sum(keep) < 5) {
      warning(sprintf("CalibrationPlot: <5 non-missing points at tau=%g.", eval_time))
      next
    }
    pred_use   <- as.numeric(prediction[keep])
    pseudo_use <- as.numeric(pseudo[keep])
    x          <- pred_use
    y          <- pseudo_use
    if (length(unique(x)) < length(x)) x <- jitter(x, factor = 1e-6)

    if (!loess_smoothing) {
      bw  <- if (is.null(bandwidth)) prodlim::neighborhood(x)$bandwidth else bandwidth
      nbh <- prodlim::meanNeighbors(x = x, y = y, bandwidth = bw)
      plotFrames[[i]] <- data.frame(pred = nbh$uniqueX, obs = nbh$averageY,
                                    lower = NA, upper = NA)
    } else {
      ord        <- order(pred_use)
      pred_use   <- pred_use[ord]
      pseudo_use <- pseudo_use[ord]
      pseu       <- data.frame(risk = pred_use, pseudovalue = pseudo_use)
      fit_loess  <- stats::loess(pseudovalue ~ risk, data = pseu,
                                  degree = 1, span = 0.3)
      sm         <- stats::predict(fit_loess, se = TRUE)
      plotFrames[[i]] <- data.frame(
        pred  = pseu$risk,
        obs   = sm$fit,
        lower = pmax(sm$fit - stats::qt(0.975, sm$df) * sm$se, 0),
        upper = pmin(sm$fit + stats::qt(0.975, sm$df) * sm$se, 1)
      )
    }

    error         <- plotFrames[[i]]$pred - plotFrames[[i]]$obs
    measures[[i]] <- data.frame(
      tau  = eval_time,
      ICI  = mean(abs(error), na.rm = TRUE),
      E50  = suppressWarnings(stats::quantile(abs(error), 0.5, na.rm = TRUE)),
      E90  = suppressWarnings(stats::quantile(abs(error), 0.9, na.rm = TRUE)),
      Emax = max(abs(error), na.rm = TRUE),
      RSB  = sqrt(mean(error^2, na.rm = TRUE))
    )

    if (graph) {
      df           <- plotFrames[[i]]
      pred_hist    <- prediction[!is.na(prediction)]
      x_max        <- max(df$pred, na.rm = TRUE)
      bin_breaks   <- seq(0, x_max, length.out = 101)
      freqs        <- table(cut(pred_hist, breaks = bin_breaks,
                                include.lowest = TRUE))
      bins         <- bin_breaks[-1]
      idx_valid    <- which(freqs > 0)
      freqs_valid  <- as.numeric(freqs[idx_valid])
      bins_valid   <- bins[idx_valid]
      spike_bounds <- c(-0.075, 0)
      spikes_df    <- NULL
      if (length(freqs_valid) > 0 && max(freqs_valid) > min(freqs_valid)) {
        fr_sc <- spike_bounds[1] +
          (spike_bounds[2] - spike_bounds[1]) *
          (freqs_valid - min(freqs_valid)) /
          (max(freqs_valid) - min(freqs_valid))
        spikes_df <- data.frame(x = bins_valid,
                                y0 = spike_bounds[1], y1 = fr_sc)
      }

      p <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = obs)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                              alpha = 0.2) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                              linetype = "dashed", colour = "red") +
        ggplot2::scale_y_continuous(breaks = seq(0, 0.6, by = 0.1),
                                    limits = c(spike_bounds[1], 0.6)) +
        ggplot2::coord_cartesian(xlim = c(0, x_max), expand = FALSE) +
        ggplot2::labs(
          x     = "Predicted risk",
          y     = "Observed proportions (pseudo-values)",
          title = paste0(model_name, " calibration for cause ", cause,
                         " at time ", round(eval_time, 1), e_txt)
        ) +
        ggplot2::theme_minimal()

      if (!is.null(spikes_df))
        p <- p + ggplot2::geom_segment(
          data        = spikes_df,
          inherit.aes = FALSE,
          ggplot2::aes(x = x, xend = x, y = y0, yend = y1)
        )
      graphs[[i]] <- p
    }
  }

  list(graphs = graphs, values = plotFrames,
       calib.measures = measures, OE_summary = OE_summary)
}


#' Observed/Expected ratio for competing risks
#'
#' @param data     Data frame with `time` and `status` columns.
#' @param pred_mat Numeric matrix `[n, Tm]` of predicted CIFs.
#' @param cause    Integer cause of interest.
#' @param times    Numeric vector of evaluation times.
#' @param horizon  Scalar evaluation time.
#'
#' @return A named numeric vector `c(OE, lower, upper)`.
#' @export
OEComputation <- function(data, pred_mat, cause, times, horizon) {
  obj <- summary(
    survival::survfit(
      survival::Surv(time, factor(status, levels = c(0L, 1L, 2L))) ~ 1,
      data = data
    ),
    times = horizon
  )
  aj          <- list(obs = obj$pstate[, cause + 1L],
                      se  = obj$std.err[, cause + 1L])
  horizon_idx <- which.min(abs(times - horizon))
  OE          <- aj$obs / mean(pred_mat[, horizon_idx])
  c(OE    = OE,
    lower = exp(log(OE) - stats::qnorm(0.975) * aj$se / aj$obs),
    upper = exp(log(OE) + stats::qnorm(0.975) * aj$se / aj$obs))
}


#' IBS wrapper for use from Python or external processes
#'
#' Computes the Integrated Brier Score for each cause in a 3-D CIF array
#' using `riskRegression::Score()`.
#'
#' @param out         3-D numeric array `[n, K, Tm]`.
#' @param e_va        Integer vector of observed event codes (length `n`).
#' @param t_va        Numeric vector of observed times (length `n`).
#' @param times       Numeric vector of evaluation times (length `Tm`).
#' @param causes      Integer vector of cause codes (length `K`).
#' @param time_col    Column name for time (default `"time"`).
#' @param status_col  Column name for status (default `"event"`).
#' @param cens.method Passed to `Score()` (default `"ipcw"`).
#' @param cens.model  Passed to `Score()` (default `"km"`).
#' @param cens.code   Censoring code (default `0`).
#' @param force_sequential If `TRUE`, force sequential execution and
#'   single-threaded BLAS (default `TRUE`).
#' @param quiet       Suppress output (default `TRUE`).
#' @param log_file    Optional path to a log file.
#' @param debug       Return debug information? (default `FALSE`).
#' @param debug_n     Number of rows to include in debug output.
#'
#' @return A numeric matrix `[K, Tm]` of IBS values, or a debug list if
#'   `debug = TRUE`.
#' @export
ibs_from_r <- function(out, e_va, t_va, times, causes,
                        time_col          = "time",
                        status_col        = "event",
                        cens.method       = "ipcw",
                        cens.model        = "km",
                        cens.code         = 0L,
                        force_sequential  = TRUE,
                        quiet             = TRUE,
                        log_file          = NULL,
                        debug             = FALSE,
                        debug_n           = 5L) {
  if (isTRUE(force_sequential)) {
    foreach::registerDoSEQ()
    data.table::setDTthreads(1L)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
  }

  log <- function(...) {
    msg <- paste0(paste(..., collapse = ""), "\n")
    if (!is.null(log_file)) cat(msg, file = log_file, append = TRUE)
    else if (!isTRUE(quiet)) cat(msg)
  }

  d  <- dim(out)
  if (length(d) != 3L) stop("`out` must be a 3D array (n, K, T).")
  n  <- d[1L]; K <- d[2L]; Tm <- d[3L]
  log("dim(out) = ", paste(d, collapse = " x "))
  if (length(times)  != Tm) stop("length(times) must equal dim(out)[3].")
  if (length(causes) != K)  stop("length(causes) must equal dim(out)[2].")

  times <- as.numeric(times)
  ord   <- order(times)
  if (any(ord != seq_along(times))) {
    out   <- out[, , ord, drop = FALSE]
    times <- times[ord]
  }
  if (any(duplicated(times))) stop("`times` must be unique.")

  test        <- data.frame(as.numeric(t_va), as.integer(e_va))
  names(test) <- c(time_col, status_col)
  f           <- stats::as.formula(
    sprintf("prodlim::Hist(%s, %s) ~ 1", time_col, status_col)
  )

  run_quiet <- function(expr) {
    if (!isTRUE(quiet)) return(eval.parent(substitute(expr)))
    res <- NULL
    invisible(utils::capture.output({
      res <- withCallingHandlers(
        eval.parent(substitute(expr)),
        warning = function(w) invokeRestart("muffleWarning")
      )
    }))
    res
  }

  ibs_mat <- matrix(NA_real_, nrow = length(causes), ncol = length(times),
                    dimnames = list(as.character(causes),
                                   as.character(times)))

  for (i in seq_len(K)) {
    k      <- as.integer(causes[i])
    M      <- out[, i, , drop = FALSE]
    dim(M) <- c(n, Tm)
    M      <- as.matrix(M)
    storage.mode(M) <- "double"
    colnames(M) <- format(times, scientific = FALSE, trim = TRUE)

    k_name     <- paste0("cif_cause_", k)
    preds_list <- stats::setNames(list(M), k_name)

    sc <- run_quiet(
      riskRegression::Score(
        object       = preds_list,
        formula      = f,
        data         = test,
        metrics      = "brier",
        summary      = "ibs",
        times        = times,
        cause        = k,
        cens.method  = cens.method,
        cens.model   = cens.model,
        cens.code    = cens.code,
        split.method = "none",
        se.fit       = FALSE,
        conf.int     = FALSE,
        progress.bar = FALSE,
        verbose      = 0,
        censoring.save.memory = TRUE,
        keep         = NULL,
        plots        = NULL
      )
    )

    ibs_val        <- sc$Brier$score[sc$Brier$score$model == k_name, "IBS"]
    ibs_mat[i, ]   <- as.numeric(ibs_val)
    rm(sc, M, preds_list); gc(FALSE)
  }

  if (isTRUE(debug)) {
    out_slice <- out[seq_len(min(n, debug_n)),
                     seq_len(min(K, debug_n)),
                     seq_len(min(Tm, debug_n)), drop = FALSE]
    return(list(
      ibs   = ibs_mat,
      debug = list(dim_out   = d, dim_slice = dim(out_slice),
                   times_head = utils::head(times, debug_n),
                   causes = causes, test_head = utils::head(test, debug_n),
                   event_table = table(test[[status_col]], useNA = "ifany"))
    ))
  }
  ibs_mat
}


#' Decision curve analysis / clinical utility plot
#'
#' Computes net benefit curves and (optionally) interventions-avoided curves
#' for one or more predicted risk columns, handling competing risks.
#'
#' @param data           Data frame with outcome columns and predictor columns.
#' @param outcome        Name of the event/status column.
#' @param ttoutcome      Name of the time column.
#' @param timepoint      Scalar evaluation time.
#' @param predictors     Character vector of predicted-risk column names
#'   (must be probabilities in [0, 1]).
#' @param xstart,xstop,xby  Threshold grid parameters.
#' @param ymin           Minimum y-axis value.
#' @param probability    Logical vector indicating which predictors are
#'   already probabilities.
#' @param harm           Numeric vector of harm values per predictor.
#' @param graph          Logical; produce a ggplot? (default `TRUE`).
#' @param intervention   Logical; plot interventions-avoided instead of
#'   net benefit? (default `FALSE`).
#' @param interventionper Denominator for the interventions-avoided scale.
#' @param smooth         Logical; apply LOESS smoothing?
#' @param loess.span     LOESS span (default 0.10).
#' @param cmprsk         Logical; use competing risks cumulative incidence
#'   (default `TRUE`).
#'
#' @return A list with `N`, `predictors`, `net.benefit`,
#'   `interventions.avoided`, `interventions.avoided.per`, and optionally
#'   `plot` (a ggplot object).
#' @export
ClinicalUtility <- function(data,
                             outcome,
                             ttoutcome,
                             timepoint,
                             predictors,
                             xstart         = 0.01,
                             xstop          = 0.99,
                             xby            = 0.01,
                             ymin           = -0.05,
                             probability    = NULL,
                             harm           = NULL,
                             graph          = TRUE,
                             intervention   = FALSE,
                             interventionper = 100,
                             smooth         = FALSE,
                             loess.span     = 0.10,
                             cmprsk         = TRUE) {

  data <- data[stats::complete.cases(
    data[c(outcome, ttoutcome, predictors)]),
    c(outcome, ttoutcome, predictors)]
  data <- as.data.frame(data)

  if (!cmprsk &&
      length(data[!(data[[outcome]] %in% c(0, 1)), outcome]) > 0)
    stop("outcome must be coded as 0 and 1 when cmprsk = FALSE.")
  if (!inherits(data, "data.frame")) stop("data must be a data.frame.")
  if (xstart < 0 || xstart > 1) stop("xstart must lie between 0 and 1.")
  if (xstop  < 0 || xstop  > 1) stop("xstop must lie between 0 and 1.")
  if (xby   <= 0 || xby    >= 1) stop("xby must lie between 0 and 1.")
  if (xstart >= xstop)            stop("xstop must be larger than xstart.")

  pred.n <- length(predictors)
  if (length(probability) > 0 && pred.n != length(probability))
    stop("Length of probability must match number of predictors.")
  if (length(harm) > 0 && pred.n != length(harm))
    stop("Length of harm must match number of predictors.")
  if (length(harm)        == 0L) harm        <- rep(0, pred.n)
  if (length(probability) == 0L) probability <- rep(TRUE, pred.n)
  if (any(predictors %in% c("all", "none")))
    stop("Predictor names cannot be 'all' or 'none'.")

  for (m in seq_len(pred.n)) {
    if (!probability[m] %in% c(TRUE, FALSE))
      stop("Each probability element must be TRUE or FALSE.")
    if (probability[m] &&
        (max(data[[predictors[m]]]) > 1 || min(data[[predictors[m]]]) < 0))
      stop(predictors[m], " must be between 0 and 1.")
    if (!probability[m]) {
      model <- survival::coxph(
        survival::Surv(data[[ttoutcome]], data[[outcome]]) ~
          data[[predictors[m]]]
      )
      surv.data <- data.frame(0)
      pred <- 1 - summary(
        survival::survfit(model, newdata = surv.data),
        time = timepoint, extend = TRUE
      )$surv
      data[[predictors[m]]] <- pred
      message(predictors[m], " converted to probability via Cox PH.")
    }
  }

  N  <- nrow(data)
  if (!cmprsk) {
    km <- survival::survfit(
      survival::Surv(data[[ttoutcome]], data[[outcome]]) ~ 1)
    pd <- 1 - summary(km, times = timepoint, extend = TRUE)$surv
  } else {
    cr <- cmprsk::cuminc(data[[ttoutcome]], data[[outcome]])
    pd <- cmprsk::timepoints(cr, times = timepoint)$est[1L]
  }

  nb    <- data.frame(threshold = seq(xstart, xstop, by = xby))
  interv <- nb
  error  <- NULL

  nb[["all"]]  <- pd - (1 - pd) * nb$threshold / (1 - nb$threshold)
  nb[["none"]] <- 0

  for (m in seq_len(pred.n)) {
    nb[[predictors[m]]] <- NA_real_
    for (t in seq_len(nrow(nb))) {
      px <- sum(data[[predictors[m]]] > nb$threshold[t]) / N
      if (px == 0) {
        error <- c(error, paste0(predictors[m], ": no obs above ",
                                  nb$threshold[t] * 100, "%."))
        break
      }
      if (!cmprsk) {
        km2 <- survival::survfit(
          survival::Surv(
            data[data[[predictors[m]]] > nb$threshold[t], ttoutcome],
            data[data[[predictors[m]]] > nb$threshold[t], outcome]
          ) ~ 1
        )
        pdg <- 1 - summary(km2, times = timepoint, extend = TRUE)$surv
        if (!length(pdg)) { error <- c(error, "no followup"); break }
      } else {
        idx <- which(data[[predictors[m]]] > nb$threshold[t])
        if (!length(idx)) { error <- c(error, "no obs"); break }
        if (!any(data[[outcome]][idx] != 0L, na.rm = TRUE)) {
          nb[t, predictors[m]] <- NA_real_; next
        }
        cr2 <- cmprsk::cuminc(data[[ttoutcome]][idx], data[[outcome]][idx])
        pdg <- cmprsk::timepoints(cr2, times = timepoint)$est[1L]
        if (is.na(pdg)) { error <- c(error, "no followup"); break }
      }
      nb[t, predictors[m]] <- pdg * px -
        (1 - pdg) * px * nb$threshold[t] / (1 - nb$threshold[t]) - harm[m]
    }
    interv[[predictors[m]]] <- (nb[[predictors[m]]] - nb[["all"]]) *
      interventionper / (interv$threshold / (1 - interv$threshold))
  }

  if (length(error))
    message(paste(error, collapse = "\n"), " Net benefit not calculable.")

  if (smooth) {
    for (m in seq_len(pred.n)) {
      valid <- !is.na(nb[[predictors[m]]])
      if (sum(valid) >= 3L) {
        lws <- stats::loess(nb[[predictors[m]]][valid] ~ nb$threshold[valid],
                            span = loess.span)
        nb[valid, paste0(predictors[m], "_sm")] <- lws$fitted
        lws <- stats::loess(interv[[predictors[m]]][valid] ~
                              interv$threshold[valid],
                            span = loess.span)
        interv[valid, paste0(predictors[m], "_sm")] <- lws$fitted
      }
    }
  }

  results <- list(
    N                    = N,
    predictors           = data.frame(predictor = predictors,
                                      harm.applied = harm,
                                      probability = probability,
                                      stringsAsFactors = FALSE),
    interventions.avoided.per = interventionper,
    net.benefit          = nb,
    interventions.avoided = interv
  )

  if (graph) {
    plot_data <- if (intervention) interv else nb
    ylab_str  <- if (intervention)
      paste("Net reduction in interventions per", interventionper, "patients")
    else
      "Net benefit"
    results$plot <- plot_dca_gg(nb        = plot_data,
                                 predictors = predictors,
                                 ylab       = ylab_str,
                                 ymin       = ymin)
  }

  results
}


#' Internal ggplot2 helper for DCA / clinical utility plots
#'
#' @param nb          Data frame with threshold + model columns (as produced
#'   by [ClinicalUtility()]).
#' @param predictors  Character vector of predictor column names.
#' @param ylab        Y-axis label.
#' @param ymin        Minimum y-axis value (passed to `coord_cartesian`).
#' @param model_cols  Named list of hex colour codes keyed by model name.
#'
#' @return A ggplot object.
#' @export
plot_dca_gg <- function(nb, predictors,
                         ylab = "Net benefit",
                         ymin = NULL,
                         model_cols = list(
                           csCPH   = "#F8766D",
                           DeepHit = "#B79F00",
                           DeSurv  = "#00BA38",
                           FGR     = "#00BFC4",
                           FGRP    = "#619CFF",
                           RSF     = "#F564E3"
                         )) {
  predictors   <- as.character(predictors)
  extract_model <- function(s) {
    m <- regmatches(s, regexpr("(csCPH|DeepHit|DeSurv|FGRP|FGR|RSF)", s))
    ifelse(nzchar(m), m, s)
  }
  label_map   <- stats::setNames(extract_model(predictors), predictors)
  pred_labels <- unname(label_map[predictors])
  model_cols  <- unlist(model_cols)
  col_map     <- c("Treat none" = "black", "Treat all" = "grey40",
                   model_cols[pred_labels])

  df <- nb |>
    dplyr::select(threshold, dplyr::all_of(c("all", "none", predictors))) |>
    tidyr::pivot_longer(-threshold, names_to = "model_raw", values_to = "y") |>
    dplyr::mutate(
      model = dplyr::case_when(
        model_raw == "none" ~ "Treat none",
        model_raw == "all"  ~ "Treat all",
        TRUE ~ dplyr::recode(model_raw, !!!label_map)
      ),
      model = factor(model, levels = c("Treat none", "Treat all", pred_labels))
    )

  p <- ggplot2::ggplot(df, ggplot2::aes(threshold, y, color = model)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position  = "right",
                   legend.title     = ggplot2::element_blank(),
                   legend.text      = ggplot2::element_text(size = 8)) +
    ggplot2::labs(x = "Threshold probability", y = ylab)

  if (!is.null(ymin))
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, NA))
  p
}

#' Censoring probabilities at subject-specific times (IPCW, individual)
#'
#' Estimates the censoring distribution via a reversed Kaplan-Meier and
#' returns subject-specific IPCW weights evaluated at each individual's
#' observed time.
#'
#' @param time Numeric vector of observed times.
#' @param status Integer vector of event/status codes.
#' @param cens.code Integer code denoting censoring (typically 0).
#' @param predictions Numeric vector of predicted risks (attached as a
#'   column in the returned data frame for convenience).
#'
#' @return A data frame with columns `time`, `status`, `predictions`,
#'   `censor.status`, and `ipcw.subject.times`.
#' @noRd
censor.prob.KM.individual <- function(time, status, cens.code, predictions) {
  tmp <- data.frame(time = time, status = status, predictions = predictions)
  tmp$censor.status <- ifelse(status == cens.code, 0, 1)
  rownames(tmp) <- seq_len(nrow(tmp))
  tmp <- tmp[order(tmp$time), ]
  fit <- prodlim::prodlim(
    formula = survival::Surv(time, censor.status) ~ 1,
    data    = tmp,
    reverse = TRUE
  )
  ipcw.subject.times <- prodlim::predictSurvIndividual(fit, lag = 0)
  tmp$ipcw.subject.times <- ipcw.subject.times
  tmp
}
