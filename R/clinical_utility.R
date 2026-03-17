#' Decision curve analysis for a single cause
#'
#' Computes net benefit curves for a single predicted risk vector (one model,
#' one cause) against an observed competing-risks outcome.  Called exclusively
#' from [compute_metrics()] when `"clinical_utility"` is in `metrics`.
#'
#' @param predictions  Numeric vector of predicted CIF values in `[0, 1]`
#'   for one cause (length `n`).
#' @param cr           A [cr_data()] object.
#' @param cause        Integer cause code.
#' @param cause_nm     Character cause name (e.g. `"cause_1"`).
#' @param timepoint    Scalar evaluation time.
#' @param xstart,xstop,xby  Threshold grid parameters.
#' @param ymin         Minimum y-axis value for the plot.
#' @param harm         Numeric scalar harm value (default `0`).
#' @param intervention Logical; plot interventions-avoided instead of net
#'   benefit? (default `FALSE`).
#' @param interventionper Denominator for the interventions-avoided scale
#'   (default `100`).
#' @param smooth       Logical; apply LOESS smoothing? (default `FALSE`).
#' @param loess.span   LOESS span (default `0.10`).
#' @param graph        Logical; produce a ggplot? (default `TRUE`).
#'
#' @return A list with elements `net_benefit` (data frame with columns
#'   `threshold`, `all`, `none`, and `model`), `interventions_avoided`
#'   (same structure), and optionally `plot` (a ggplot object).
#' @noRd
.clinical_utility <- function(predictions,
                              cr,
                              cause,
                              cause_nm,
                              timepoint,
                              xstart         = 0.01,
                              xstop          = 0.99,
                              xby            = 0.01,
                              ymin           = -0.05,
                              harm           = 0,
                              intervention   = FALSE,
                              interventionper = 100,
                              smooth         = FALSE,
                              loess.span     = 0.10,
                              graph          = TRUE) {
  
  time   <- cr@data[[cr@time_var]]
  status <- cr@data[[cr@event_var]]
  N      <- length(predictions)
  
  # Prevalence at timepoint via competing-risks cumulative incidence
  cr_fit <- cmprsk::cuminc(time, status)
  pd     <- cmprsk::timepoints(cr_fit, times = timepoint)$est[1L]
  
  nb     <- data.frame(threshold = seq(xstart, xstop, by = xby))
  interv <- nb
  error  <- NULL
  
  nb[["all"]]   <- pd - (1 - pd) * nb$threshold / (1 - nb$threshold)
  nb[["none"]]  <- 0
  nb[["model"]] <- NA_real_
  
  for (t in seq_len(nrow(nb))) {
    px <- sum(predictions > nb$threshold[t]) / N
    if (px == 0) {
      error <- c(error, sprintf(
        "%s: no observations above threshold %.0f%%.", cause_nm,
        nb$threshold[t] * 100
      ))
      break
    }
    idx <- which(predictions > nb$threshold[t])
    if (!any(status[idx] != cr@cens_code, na.rm = TRUE)) {
      nb[t, "model"] <- NA_real_
      next
    }
    cr2 <- cmprsk::cuminc(time[idx], status[idx])
    pdg <- cmprsk::timepoints(cr2, times = timepoint)$est[1L]
    if (is.na(pdg)) { nb[t, "model"] <- NA_real_; next }
    nb[t, "model"] <- pdg * px -
      (1 - pdg) * px * nb$threshold[t] / (1 - nb$threshold[t]) - harm
  }
  
  interv[["all"]]   <- nb[["all"]]
  interv[["none"]]  <- 0
  interv[["model"]] <- (nb[["model"]] - nb[["all"]]) *
    interventionper / (interv$threshold / (1 - interv$threshold))
  
  if (length(error))
    message(paste(error, collapse = "\n"), " Net benefit not calculable.")
  
  if (smooth) {
    valid <- !is.na(nb[["model"]])
    if (sum(valid) >= 3L) {
      lws <- stats::loess(nb[["model"]][valid] ~ nb$threshold[valid],
                          span = loess.span)
      nb[valid,     "model_sm"] <- lws$fitted
      lws <- stats::loess(interv[["model"]][valid] ~ interv$threshold[valid],
                          span = loess.span)
      interv[valid, "model_sm"] <- lws$fitted
    }
  }
  
  results <- list(
    net_benefit           = nb,
    interventions_avoided = interv
  )
  
  if (graph) {
    plot_data <- if (intervention) interv else nb
    ylab_str  <- if (intervention)
      paste("Net reduction in interventions per", interventionper, "patients")
    else
      "Net benefit"
    results$plot <- .plot_dca_gg(
      nb        = plot_data,
      cause_nm  = cause_nm,
      timepoint = timepoint,
      ylab      = ylab_str,
      ymin      = ymin
    )
  }
  
  results
}


#' Internal ggplot2 helper for DCA / clinical utility plots
#'
#' @param nb        Data frame with columns `threshold`, `all`, `none`,
#'   `model` (as produced by [.clinical_utility()]).
#' @param cause_nm  Character cause name used in the plot title.
#' @param timepoint Scalar evaluation time used in the plot title.
#' @param ylab      Y-axis label.
#' @param ymin      Minimum y-axis value (passed to `coord_cartesian`).
#'
#' @return A ggplot object.
#' @noRd
.plot_dca_gg <- function(nb,
                         cause_nm,
                         timepoint,
                         ylab = "Net benefit",
                         ymin = NULL) {
  
  df <- nb |>
    dplyr::select(threshold, dplyr::all_of(c("all", "none", "model"))) |>
    tidyr::pivot_longer(-threshold, names_to = "curve", values_to = "y") |>
    dplyr::mutate(
      curve = dplyr::case_when(
        curve == "none"  ~ "Treat none",
        curve == "all"   ~ "Treat all",
        curve == "model" ~ cause_nm
      ),
      curve = factor(curve, levels = c("Treat none", "Treat all", cause_nm))
    )
  
  col_map <- c("Treat none" = "black", "Treat all" = "grey40",
               stats::setNames("#00BFC4", cause_nm))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(threshold, y, color = curve)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   legend.title    = ggplot2::element_blank(),
                   legend.text     = ggplot2::element_text(size = 8)) +
    ggplot2::labs(
      x     = "Threshold probability",
      y     = ylab,
      title = paste0("Clinical utility for ", cause_nm,
                     " at time ", round(timepoint, 1))
    )
  
  if (!is.null(ymin))
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, NA))
  p
}