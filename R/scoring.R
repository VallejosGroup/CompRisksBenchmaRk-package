#' @title Scoring and Results Utilities
#' @description Functions for computing, summarising, aggregating, and
#'   plotting out-of-sample benchmark results.
#' @name scoring
NULL



#' Extract penalised Fine-Gray regression coefficients
#'
#' Convenience helper to pull the named coefficient vector from a fitted
#' FGRP model object (as returned by the registered `"FGRP"` model).
#'
#' @param fit_obj A fit object returned by the FGRP `fit()`
#'   function.
#' @param tol     Threshold below which absolute coefficient values are
#'   considered zero (default `1e-8`).
#'
#' @return A data frame with columns `cause`, `term`, `beta`, and `kept`.
#' @export
extract_fgrp_coefs <- function(fit_obj, tol = 1e-8) {
  get_nm <- function(fp, p) {
    nm <- names(fp[["coef"]])
    if (is.null(nm)) nm <- fp[["coefnames"]]
    if (is.null(nm)) nm <- fp[["varnames"]]
    if (is.null(nm)) nm <- paste0("V", seq_len(p))
    nm
  }
  res <- lapply(seq_along(fit_obj$causes), function(i) {
    fp   <- fit_obj$fits[[i]]$fp
    beta <- as.vector(fp[["coef"]])
    nm   <- get_nm(fp, length(beta))
    data.frame(cause = fit_obj$causes[[i]], term = nm, beta = beta,
               kept = abs(beta) > tol, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}


#' Load and score benchmark predictions across models
#'
#' Reads pre-computed CIF Parquet files from `data_root/dataset/<model>/outer_v/`,
#' scores them, plots calibration curves, and returns all results in a single list.
#'
#' @param data_root  Root directory containing dataset subdirectories.
#' @param dataset    Name of the dataset subdirectory.
#' @param dataset_prior_imputation Name of the directory holding the original
#'   (pre-imputation) splits used to locate the test sets.  Required when
#'   `dataset` contains `"imputed"`.
#' @param python_models Character vector of model keys whose CIF arrays are
#'   stored in `[K, Tm, n]` order (Python convention).
#' @param r_models  Character vector of model keys whose CIF arrays are in
#'   `[n, K, Tm]` order.
#' @param outer_folds Number of outer folds.
#' @param horizon   Scalar evaluation time for calibration and O/E.
#' @param seed      Random seed for sampling example curves.
#' @param cause_of_interest Integer cause index.
#' @param sample_curves If `TRUE`, print example CIF trajectories per fold.
#'
#' @return A list with `all_test_folds`, `predictor_cols`, `all_cifs`,
#'   `all_rmlt`, `res`, `res_rmlt`, and `oes`.
#' @export
get_results <- function(data_root,
                        dataset,
                        dataset_prior_imputation = NULL,
                        python_models            = NULL,
                        r_models,
                        outer_folds,
                        horizon,
                        seed,
                        cause_of_interest,
                        sample_curves) {
  res        <- list()
  res_rmlt   <- list()
  all_cifs   <- list()
  all_test   <- list()
  all_rmlt   <- list()
  oes        <- list()
  all_test_folds <- NULL

  is_imputed <- grepl("imputed", dataset, fixed = TRUE)
  if (is_imputed && is.null(dataset_prior_imputation))
    stop("dataset_prior_imputation must be provided.", call. = FALSE)

  models <- c(python_models, r_models)

  for (model in models) {
    fold_dir <- file.path(data_root, dataset, model)
    res[[model]]      <- vector("list", outer_folds)
    res_rmlt[[model]] <- vector("list", outer_folds)

    for (v in seq_len(outer_folds)) {
      cif_path  <- file.path(fold_dir, sprintf("outer_%d", v), "cif.parquet")
      test_path <- if (!is.null(dataset_prior_imputation))
        file.path(data_root, dataset_prior_imputation,
                  sprintf("outer_%d", v), "test.parquet")
      else
        file.path(data_root, dataset,
                  sprintf("outer_%d", v), "test.parquet")

      df_cif  <- arrow::read_parquet(cif_path)
      df_test <- arrow::read_parquet(test_path)

      times     <- sort(unique(df_cif$time))
      causes    <- sort(unique(df_cif$cause))
      n         <- nrow(df_test)
      K         <- length(causes)
      Tm        <- length(times)
      cause_idx <- cause_of_interest

      if (!("row_id" %in% names(df_cif))) {
        stopifnot(n * K * Tm == nrow(df_cif))
        df_cif$row_id <- rep(df_test$row_id, each = Tm * K)
      }

      if (model %in% python_models) {
        ord    <- order(match(df_cif$row_id, df_test$row_id),
                        df_cif$time, df_cif$cause)
        df_cif <- df_cif[ord, , drop = FALSE]
        stopifnot(length(df_cif$cif) == n * K * Tm)
        out <- aperm(array(df_cif$cif, dim = c(K, Tm, n)), c(3, 1, 2))
      } else {
        out <- array(df_cif$cif, dim = c(n, K, Tm))
      }

      if (sample_curves) {
        set.seed(seed)
        ids <- sample.int(n, min(20, n))
        mat <- out[ids, cause_idx, ]
        p   <- ggplot2::ggplot(
          data.frame(
            time = rep(times, each = length(ids)),
            id   = rep(ids,   times = length(times)),
            cif  = as.vector(mat)
          ),
          ggplot2::aes(time, cif, group = id)
        ) +
          ggplot2::geom_line(color = "steelblue", alpha = 0.35,
                             linewidth = 0.4) +
          ggplot2::coord_cartesian(ylim = c(0, 1)) +
          ggplot2::labs(
            x     = "Time",
            y     = paste0("CIF (cause ", causes[cause_idx], ")"),
            title = paste("Fold", v, "- Cause", causes[cause_idx],
                          "- Model", model)
          ) +
          ggplot2::theme_minimal(base_size = 11)
        print(p)
      }

      cr_test <- cr_data(df_test, time_var = "time", event_var = "event")

      cif_out <- list(cif = out, time_grid = times)

      res[[model]][[v]] <- compute_metrics(
        cr             = cr_test,
        eval_times     = times,
        cif            = cif_out,
        metrics        = c("Brier", "IBS", "tdAUC", "calib_measures"),
        args_riskRegression = list(cens.method = "ipcw", cens.model = "km"),
        collapse_as_df = FALSE
      )

      all_cifs[[model]][[v]] <- out
      all_test[[v]]          <- df_test

      all_rmlt[[model]][[v]] <- compute_rmlt(cr_test,
                                             cif  = cif_out,
                                             maxT = max(times))
      res_rmlt[[model]][[v]] <- compute_metrics(
        cr             = cr_test,
        eval_times     = max(times),
        cif            = cif_out,
        metrics        = "cindex_rmlt",
        collapse_as_df = FALSE
      )
    }

    all_cifs[[model]] <- abind::abind(all_cifs[[model]], along = 1)
    if (is.null(all_test_folds))
      all_test_folds <- do.call(rbind, all_test)

    cause1      <- all_cifs[[model]][, cause_of_interest, ]
    horizon_idx <- which.min(abs(times - horizon))
    horizon_use <- times[horizon_idx]

    plotFrame <- CalibrationPlot(
      model_name       = model,
      predictions      = cause1[, horizon_idx],
      data             = all_test_folds,
      time             = all_test_folds$time,
      status           = all_test_folds$event,
      tau              = horizon_use,
      cause            = cause_of_interest,
      cens.code        = 0,
      loess_smoothing  = TRUE,
      predictions.type = "CIF",
      graph            = TRUE
    )
    print(plotFrame$graphs)

    col_name <- paste0("cif_cause", cause_of_interest, "_", model,
                       "_t", round(horizon_use, 1))
    all_test_folds[[col_name]] <- as.numeric(cause1[, horizon_idx])
    oes[[col_name]]            <- plotFrame$OE_summary
  }

  predictor_cols <- paste0("cif_cause", cause_of_interest, "_",
                           models, "_t", round(horizon, 1))

  list(all_test_folds = all_test_folds,
       predictor_cols = predictor_cols,
       all_cifs       = all_cifs,
       all_rmlt       = all_rmlt,
       res            = res,
       res_rmlt       = res_rmlt,
       oes            = oes)
}


#' Summarise nested cross-validation results across outer folds
#'
#' @param results A list of per-fold result lists as returned by
#'   [nested_cv_from_bench()].
#' @param num_causes Integer; number of competing causes.
#' @param times Numeric vector of evaluation time points.
#'
#' @return A data frame with one row per model × fold × cause × time,
#'   containing Brier score, IBS, AUC, C-index (pec),
#'   and calibration measures.
#' @export
summarize_out_of_sample <- function(results, num_causes, times) {
  tab <- list()

  for (i in seq_along(results)) {
    r      <- results[[i]]
    model  <- names(results)[[i]]
    outer_folds <- length(r)

    for (v in seq_len(outer_folds)) {
      res_fold <- r[[v]]
      causes   <- num_causes

      get_or_na <- function(x, len)
        if (!is.null(x)) x else array(NA_real_, len)

      for (cause in seq_len(causes)) {
        bs         <- get_or_na(res_fold$Brier[[cause]],         length(times))
        ibs        <- get_or_na(res_fold$IBS[[cause]],           length(times))
        auc        <- get_or_na(res_fold$tdAUC[[cause]],         length(times))
        cidx_pec   <- get_or_na(res_fold$cindex_t_year[[cause]], length(times))
        calib_meas <- if (!is.null(res_fold$calib_measures[[cause]]))
          res_fold$calib_measures[[cause]]
        else
          list(ICI = NA_real_, E50 = NA_real_, E90 = NA_real_,
               Emax = NA_real_, RSB = NA_real_)

        tab[[length(tab) + 1]] <- data.frame(
          model      = model,
          fold       = v,
          cause      = cause,
          bs         = bs,
          ibs        = ibs,
          auc        = auc,
          cidx_pec   = cidx_pec,
          ICI        = calib_meas[["ICI"]],
          E50        = calib_meas[["E50"]],
          E90        = calib_meas[["E90"]],
          Emax       = calib_meas[["Emax"]],
          RSB        = calib_meas[["RSB"]],
          times      = times
        )
      }
    }
  }
  do.call(rbind, tab)
}


#' Aggregate out-of-sample metrics across folds
#'
#' For each metric, computes mean, median, SD, and 95% quantile interval
#' over folds, grouped by model, cause, and time.
#'
#' @param df      Data frame as returned by [summarize_out_of_sample()].
#' @param metrics Character vector of column names to aggregate.
#'
#' @return A named list (one element per metric) of aggregated data frames.
#' @export
aggregate_out_of_sample <- function(df, metrics) {
  agg_one <- function(metric) {
    aggregate(
      df[[metric]],
      by  = list(model = df$model, cause = df$cause, time = df$time),
      FUN = function(x)
        c(mean   = mean(x, na.rm = TRUE),
          median = stats::median(x, na.rm = TRUE),
          sd     = stats::sd(x, na.rm = TRUE),
          n      = sum(!is.na(x)),
          lwr    = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
          upr    = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE))
    )
  }
  out <- lapply(metrics, function(m) data.frame(agg_one(m)))
  names(out) <- metrics
  out
}


#' Plot aggregated out-of-sample metrics
#'
#' Produces one ggplot per requested metric showing median ± 95%
#' interval trajectories over time, coloured by model.
#'
#' @param agg     Output of [aggregate_out_of_sample()].
#' @param metrics Character vector of metrics to plot.
#' @param cause   Integer cause index to filter on (default: all causes).
#' @param plot_tau If `TRUE`, adds reference lines at evaluation time points.
#'
#' @return Called for its side effect (prints plots).
#' @export
plot_out_of_sample <- function(agg,
                                metrics  = c("bs", "ibs", "auc",
                                             "cidx_pec",
                                             "ICI", "E50", "E90",
                                             "Emax", "RSB"),
                                cause    = NULL,
                                plot_tau = TRUE) {
  for (metric in metrics) {
    ylimit <- c(0, 1)
    df     <- agg[[metric]]
    if (!is.null(cause)) df <- df[df$cause == cause, , drop = FALSE]

    tau <- NULL

    if (metric %in% c("bs", "auc", "cidx_pec",
                       "ICI", "E50", "E90", "Emax", "RSB")) {
      tau <- df$time
      if (length(unique(tau)) == 1 && metric == "cidx_pec") {
        df$model <- factor(df$model, levels = unique(df$model))
        p <- ggplot2::ggplot(df,
                             ggplot2::aes(x = model,
                                         y = .data[["x"]][, "median"],
                                         color = model)) +
          ggplot2::geom_point(size = 3) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data[["x"]][, "lwr"],
                         ymax = .data[["x"]][, "upr"]),
            width = 0.15, linewidth = 1
          ) +
          ggplot2::ylim(ylimit) +
          ggplot2::labs(x = "Model", y = toupper(metric), color = "Model") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(legend.position = "top")
        print(p)
        next
      }
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time,
                                          y = .data[["x"]][, "median"],
                                          group = model)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["x"]][, "lwr"],
                     ymax = .data[["x"]][, "upr"],
                     fill = model),
        alpha = 0.15, color = NA
      ) +
      ggplot2::geom_line(ggplot2::aes(color = model), linewidth = 0.9) +
      ggplot2::ylim(ylimit) +
      ggplot2::labs(x = "Time", y = toupper(metric),
                    color = "Model", fill = "Model") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "top")

    if (plot_tau && metric %in% c("bs", "auc", "cidx_pec",
                                   "ICI", "E50", "E90", "Emax", "RSB") &&
        !is.null(tau)) {
      p <- p + ggplot2::geom_vline(xintercept = tau,
                                    linetype = "dashed", color = "grey70")
    }
    print(p)
  }
}


#' Summarise outer results from [nested_cv_from_bench()]
#'
#' Extracts best hyperparameter configurations and IBS trajectories from
#' the per-fold result list.
#'
#' @param results List returned by [nested_cv_from_bench()].
#'
#' @return A list with elements `configs_per_cause`, `ibs_per_fold_cause`,
#'   `ibs_per_time`, and `ibs_summary`.
#' @export
summarize_outer_results <- function(results) {
  cfg_rows      <- list()
  ibs_rows      <- list()
  ibs_time_rows <- list()

  for (i in seq_along(results)) {
    r      <- results[[i]]
    fold   <- r$fold
    model  <- r$model_key
    times  <- as.numeric(r$times)
    causes <- r$metrics$causes

    for (ci in seq_along(causes)) {
      k  <- causes[ci]
      nm <- paste0("cause_", k)

      cfg_i  <- r$best_cfg[[nm]]
      cfg_df <- cbind(
        data.frame(fold = fold, model = model, cause = k,
                   stringsAsFactors = FALSE),
        as.data.frame(as.list(cfg_i), stringsAsFactors = FALSE)
      )
      cfg_rows[[length(cfg_rows) + 1]] <- cfg_df

      ibs_vec <- as.numeric(r$metrics$ibs_scores[[ci]])
      ibs_rows[[length(ibs_rows) + 1]] <- data.frame(
        fold = fold, model = model, cause = k,
        ibs  = utils::tail(ibs_vec, 1),
        stringsAsFactors = FALSE
      )

      Tm     <- length(ibs_vec)
      t_used <- if (length(times) == Tm) times else seq_len(Tm)
      ibs_time_rows[[length(ibs_time_rows) + 1]] <- data.frame(
        fold = fold, model = model, cause = k,
        time = t_used, ibs = ibs_vec,
        stringsAsFactors = FALSE
      )
    }
  }

  rbindlist <- function(lst)
    do.call(rbind, lapply(lst, as.data.frame))

  configs_per_cause  <- rbindlist(cfg_rows)
  ibs_per_fold_cause <- rbindlist(ibs_rows)
  ibs_per_time       <- rbindlist(ibs_time_rows)

  ibs_summary <- if (nrow(ibs_per_fold_cause) > 0) {
    key <- interaction(ibs_per_fold_cause$model,
                       ibs_per_fold_cause$cause, drop = TRUE)
    spl <- split(ibs_per_fold_cause$ibs, key)
    stats <- lapply(spl, function(x)
      c(mean = mean(x), sd = stats::sd(x), median = stats::median(x)))
    out <- do.call(rbind, stats)
    keys  <- names(spl)
    parts <- read.table(text = keys, sep = ".",
                         stringsAsFactors = FALSE,
                         col.names = c("model", "cause"))
    cbind(parts, as.data.frame(out, row.names = NULL))
  } else {
    data.frame(model = character(), cause = integer(),
               mean = numeric(), sd = numeric(), median = numeric())
  }

  list(configs_per_cause  = configs_per_cause,
       ibs_per_fold_cause = ibs_per_fold_cause,
       ibs_per_time       = ibs_per_time,
       ibs_summary        = ibs_summary)
}
