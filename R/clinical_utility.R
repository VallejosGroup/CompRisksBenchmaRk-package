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
clinical_utility <- function(data,
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
#'   by [clinical_utility()]).
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
