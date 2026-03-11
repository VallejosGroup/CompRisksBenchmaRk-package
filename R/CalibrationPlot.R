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
      
      if (!all(is.na(df$lower)))
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2
        )
      
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