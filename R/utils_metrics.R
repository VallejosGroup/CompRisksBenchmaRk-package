#' @title Internal metric utilities
#' @description Internal helper functions shared across metric computations.
#' @name utils_metrics
#' @keywords internal
NULL


#' Check that `cr` is a cr_data object
#'
#' @param cr Object to check.
#' @noRd
.check_cr <- function(cr) {
  if (!methods::is(cr, "cr_data"))
    stop("`cr` must be a cr_data object.", call. = FALSE)
  invisible(NULL)
}


#' Resolve and validate the `cif` / `fit` / `cif_time_grid` arguments
#'
#' Enforces mutual exclusion of `cif` and `fit`, checks that `cif_time_grid`
#' is `NULL` when `cif` is supplied directly, and calls [predict_cif()] when
#' `fit` is supplied.  Returns the resolved `cif` list.
#'
#' @param cif           A predict_cif() list or `NULL`.
#' @param fit           A fitted model object or `NULL`.
#' @param cif_time_grid Numeric time grid (required when `fit` is non-`NULL`).
#' @param cr            A `cr_data` object (passed to [predict_cif()]).
#'
#' @return The resolved `cif` list (as returned by [predict_cif()]).
#' @noRd
.resolve_cif <- function(cif, fit, cif_time_grid, cr) {
  if (is.null(cif) && is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(fit))
    stop("Exactly one of `cif` or `fit` must be non-NULL; both are non-NULL.",
         call. = FALSE)
  if (!is.null(cif) && !is.null(cif_time_grid))
    stop("`cif_time_grid` must be NULL when `cif` is supplied directly.",
         call. = FALSE)
  if (!is.null(fit)) {
    if (is.null(cif_time_grid))
      stop("`cif_time_grid` must be provided when `fit` is non-NULL.",
           call. = FALSE)
    cif <- predict_cif(fit, newdata = cr, time_grid = cif_time_grid)
  }
  cif
}


#' Validate a cif list and unpack into array + time grid
#'
#' Checks that `cif` is a list with a 3-D numeric `$cif` array and a
#' non-empty numeric `$time_grid`, and that the number of cause dimensions
#' matches `cr@causes`.  Returns a list with elements `cif` (the array) and
#' `time_grid` (the numeric vector).
#'
#' @param cif A list as returned by [predict_cif()].
#' @param cr  A `cr_data` object.
#'
#' @return A list with elements `cif` (3-D array), `time_grid` (numeric
#'   vector), and `model_key` (character, `NULL` if not present in input).
#' @noRd
.validate_and_unpack_cif <- function(cif, cr) {
  if (!is.list(cif) || !all(c("cif", "time_grid") %in% names(cif)))
    stop("`cif` must be a list with elements `$cif` and `$time_grid`, ",
         "as returned by predict_cif().", call. = FALSE)
  if (!is.numeric(cif$cif) || length(dim(cif$cif)) != 3)
    stop("`cif$cif` must be a 3-D numeric array `[n, K, Tm]`.", call. = FALSE)
  if (!is.numeric(cif$time_grid) || length(cif$time_grid) == 0)
    stop("`cif$time_grid` must be a non-empty numeric vector.", call. = FALSE)
  if (dim(cif$cif)[2] != length(cr@causes))
    stop(sprintf(
      "`cif` has %d cause dimension(s) but `cr` contains %d cause(s) (%s).",
      dim(cif$cif)[2], length(cr@causes), paste(cr@causes, collapse = ", ")
    ), call. = FALSE)
  list(cif = cif$cif, time_grid = cif$time_grid, model_key = cif$model_key)
}


#' Trapezoidal numerical integration
#'
#' @param x Numeric vector of x-values (must be sorted).
#' @param y Numeric vector of y-values (same length as `x`).
#'
#' @return A single numeric value.
#' @noRd
.trapezoidal_integration <- function(x, y) {
  if (length(x) != length(y))
    stop("`x` and `y` must have the same length.")
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}


#' Calibration for a single cause at all tau
#'
#' @param pred_at_tau  Numeric matrix `[n, length(tau)]` of predicted CIF
#'   values for this cause, already snapped to `tau`.
#' @param tau          Numeric vector of evaluation times.
#' @param cause        Integer cause code.
#' @param cause_idx    Integer position of this cause in `cr@causes`.
#' @param cause_nm     Character name of this cause (e.g. `"cause_1"`).
#' @param model_name   Character label used in plot titles.
#' @param time         Numeric vector of observed times.
#' @param status       Integer vector of event codes.
#' @param cens_code    Integer censoring code.
#' @param cr           A [cr_data()] object.
#' @param margFit      A pre-fitted `prodlim` marginal model (loop-invariant).
#' @param loess_smoothing Logical.
#' @param bandwidth    Optional numeric bandwidth.
#' @param graph        Logical; produce ggplot objects?
#'
#' @return A list with elements `graphs`, `values`, `calib_measures`,
#'   `OE_summary` — all for this cause only.
#' @noRd
.calibration_one_cause <- function(pred_at_tau, tau, cause, cause_idx,
                                   cause_nm, model_name,
                                   time, status, cens_code, cr,
                                   margFit, loess_smoothing, bandwidth,
                                   graph) {
  # O/E at the last (or only) horizon
  horizon <- tau[length(tau)]
  obj     <- summary(
    survival::survfit(
      survival::Surv(time, factor(status, levels = c(cens_code, cr@causes))) ~ 1
    ),
    times = horizon
  )
  aj          <- list(obs = as.numeric(obj$pstate[, cause_idx + 1L]),
                      se  = as.numeric(obj$std.err[, cause_idx + 1L]))
  horizon_idx <- which.min(abs(tau - horizon))
  OE          <- aj$obs / mean(pred_at_tau[, horizon_idx])
  OE_summary  <- data.frame(
    cause = cause_nm,
    OE    = OE,
    lower = exp(log(OE) - stats::qnorm(0.975) * aj$se / aj$obs),
    upper = exp(log(OE) + stats::qnorm(0.975) * aj$se / aj$obs),
    row.names = NULL
  )
  e_txt <- sprintf(", OE = %.3f (CI: %.3f\u2013%.3f)",
                   OE_summary$OE, OE_summary$lower, OE_summary$upper)
  
  plotFrames <- vector("list", length(tau))
  measures   <- vector("list", length(tau))
  graphs     <- vector("list", length(tau))
  
  for (i in seq_along(tau)) {
    prediction <- pred_at_tau[, i]
    eval_time  <- tau[i]
    
    if (length(unique(stats::na.omit(prediction))) <= 1) {
      measures[[i]] <- data.frame(cause = cause_nm, tau = eval_time,
                                  ICI = NA, E50 = NA, E90 = NA,
                                  Emax = NA, RSB = NA)
      next
    }
    
    pseudo <- prodlim::jackknife(margFit, cause = cause, times = eval_time)
    
    keep <- !is.na(prediction) & !is.na(pseudo)
    if (sum(keep) < 5) {
      warning(sprintf(
        ".compute_calibration: <5 non-missing points at tau=%g for cause %s.",
        eval_time, cause_nm
      ))
      next
    }
    pred_use   <- as.numeric(prediction[keep])
    pseudo_use <- as.numeric(pseudo[keep])
    x          <- pred_use
    y          <- pseudo_use
    if (length(unique(x)) < length(x)) x <- jitter(x, factor = 1e-6)
    
    if (!loess_smoothing) {
      bw  <- if (is.null(bandwidth)) prodlim::neighborhood(x)$bandwidth
      else bandwidth
      nbh <- prodlim::meanNeighbors(x = x, y = y, bandwidth = bw)
      plotFrames[[i]] <- data.frame(pred = nbh$uniqueX, obs = nbh$averageY,
                                    lower = NA_real_, upper = NA_real_)
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
      cause = cause_nm,
      tau   = eval_time,
      ICI   = mean(abs(error), na.rm = TRUE),
      E50   = suppressWarnings(stats::quantile(abs(error), 0.5, na.rm = TRUE)),
      E90   = suppressWarnings(stats::quantile(abs(error), 0.9, na.rm = TRUE)),
      Emax  = max(abs(error), na.rm = TRUE),
      RSB   = sqrt(mean(error^2, na.rm = TRUE)),
      row.names = NULL
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
  
  list(
    graphs         = graphs,
    values         = plotFrames,
    calib_measures = do.call(rbind, Filter(Negate(is.null), measures)),
    OE_summary     = OE_summary
  )
}


#' Compute calibration plots and measures for all competing causes
#'
#' Internal function producing pseudo-value calibration plots and numerical
#' calibration measures for all causes in a [cr_data()] object.
#'
#' @param cr            A [cr_data()] object.
#' @param cif           A list as returned by [predict_cif()].
#'   Mutually exclusive with `fit`.
#' @param fit           A fitted model object. Mutually exclusive with `cif`.
#' @param cif_time_grid Numeric time grid; required when `fit` is supplied,
#'   must be `NULL` when `cif` is supplied.
#' @param horizon       Numeric scalar or vector of evaluation times (sorted,
#'   no NAs).  Each value is mapped to the nearest point on the CIF time grid.
#' @param tol           Maximum allowable absolute difference between a
#'   requested `horizon` value and its nearest point on the CIF time grid.
#'   An error is raised if any mapping exceeds this tolerance (default `Inf`,
#'   i.e. no check).
#' @param loess_smoothing Logical; use LOESS smoothing (default `TRUE`).
#' @param bandwidth     Optional bandwidth for nearest-neighbour smoothing.
#' @param graph         Logical; produce ggplot objects? (default `TRUE`).
#'
#' @return A list with elements `graphs`, `values`, `calib_measures`, and
#'   `OE_summary`.
#' @noRd
.compute_calibration <- function(cr,
                                 cif           = NULL,
                                 fit           = NULL,
                                 cif_time_grid = NULL,
                                 horizon,
                                 tol             = Inf,
                                 loess_smoothing = TRUE,
                                 bandwidth       = NULL,
                                 graph           = TRUE) {
  .check_cr(cr)
  
  if (missing(horizon) || !is.numeric(horizon) || length(horizon) == 0)
    stop("`horizon` must be a non-empty numeric vector.", call. = FALSE)
  if (anyNA(horizon) || is.unsorted(horizon))
    stop("`horizon` must be sorted and contain no NA values.", call. = FALSE)
  
  cif_in     <- .resolve_cif(cif, fit, cif_time_grid, cr)
  unpacked   <- .validate_and_unpack_cif(cif_in, cr)
  cif_arr    <- unpacked$cif
  model_name <- unpacked$model_key
  
  # Map each horizon value to the nearest point on the CIF time grid.
  # An error is raised if any difference exceeds `tol`.
  time_grid          <- unpacked$time_grid
  map_horizon_to_grid <- vapply(
    horizon,
    function(h) which.min(abs(time_grid - h)),
    integer(1)
  )
  diffs <- abs(horizon - time_grid[map_horizon_to_grid])
  if (any(diffs > tol))
    stop(sprintf(
      "horizon value(s) %s are more than tol=%g away from the nearest ",
      "CIF time grid point (max diff = %g). ",
      "Adjust `horizon` or increase `tol`.",
      paste(horizon[diffs > tol], collapse = ", "), tol, max(diffs)
    ), call. = FALSE)
  
  causes    <- cr@causes
  cause_nms <- paste0("cause_", causes)
  time      <- cr@data[[cr@time_var]]
  status    <- cr@data[[cr@event_var]]
  cens_code <- cr@cens_code
  
  # Fit marginal model once — shared across all causes and all horizon values
  margForm <- prodlim::Hist(time, status, cens.code = cens_code) ~ 1
  margFit  <- prodlim::prodlim(margForm, data = cr@data)
  
  graphs_out          <- stats::setNames(vector("list", length(causes)), cause_nms)
  values_out          <- stats::setNames(vector("list", length(causes)), cause_nms)
  calib_measures_list <- vector("list", length(causes))
  OE_list             <- vector("list", length(causes))
  
  for (i in seq_along(causes)) {
    k        <- causes[i]
    cause_nm <- cause_nms[i]
    
    # Extract predicted CIF for this cause at the mapped horizon indices.
    # cif_arr is [n, K, Tm]; indexing the K and Tm dims gives [n, 1, length(horizon)],
    # which is then squeezed to [n, length(horizon)].
    cif_at_horizon <- cif_arr[, i, map_horizon_to_grid, drop = FALSE][, 1, , drop = FALSE]
    dim(cif_at_horizon) <- c(dim(cif_arr)[1], length(horizon))
    
    res <- .calibration_one_cause(
      pred_at_tau     = cif_at_horizon,
      tau             = horizon,
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