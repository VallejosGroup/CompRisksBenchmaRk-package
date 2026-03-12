#' @title Performance Metrics for Competing Risks Models
#' @description Functions for computing discrimination, calibration, and
#'   clinical utility metrics for competing risks survival models.
#' @name metrics
NULL


#' Observed/Expected ratio
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