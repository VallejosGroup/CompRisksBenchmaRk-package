#' @title Nested Cross-Validation for Competing Risks Models
#' @description The core engine for running nested cross-validation over a
#'   pre-built benchmark directory created by [create_nested_splits()].
#' @name nested_cv
NULL



# ---- main function -----------------------------------------------------------

#' Run nested cross-validation from a benchmark directory
#'
#' Reads pre-built nested CV splits from disk (as created by
#' [create_nested_splits()]), fits the requested model in each outer fold
#' (with optional inner-loop hyperparameter tuning), and returns per-fold
#' performance metrics.
#'
#' For models with `needs_tuning = TRUE`, the inner loop evaluates each
#' hyperparameter configuration in `grid` using the Integrated Brier Score
#' and selects the best configuration per cause before refitting on the full
#' outer training set.
#'
#' @param out_dir   Path to the benchmark directory (default
#'   `"../BenchResults"`).
#' @param model_key Character key of a model registered with
#'   [register_cr_model()].
#' @param grid      For tunable models: a named list of hyperparameter
#'   sequences to be expanded via `tidyr::expand_grid()`.  Ignored for
#'   models with `needs_tuning = FALSE`.
#' @param seed      Base random seed (outer fold `v` uses `seed + v`).
#' @param verbose   Print progress messages.
#' @param drop_cols_after_impute Optional character vector of feature columns
#'   to exclude from the model even after they have been used in imputation.
#'
#' @return A list of length `outer_folds`.  Each element is a named list
#'   containing `fold`, `model_key`, `times`, `best_cfg` (tunable models
#'   only), and `metrics` (`causes`, `cause_bs`, `cause_ibs`).
#' @export
nested_cv_from_bench <- function(out_dir                = "../BenchResults",
                                  model_key,
                                  grid                   = NULL,
                                  seed                   = 123L,
                                  verbose                = TRUE,
                                  drop_cols_after_impute = NULL) {

  manifest_path <- file.path(out_dir, "manifest.json")
  times_path    <- file.path(out_dir, "times.json")
  meta_path     <- file.path(out_dir, "cr_metadata.json")

  meta   <- jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
  times  <- sort(unique(as.numeric(
    jsonlite::fromJSON(times_path, simplifyVector = TRUE)
  )))
  cr_meta <- if (file.exists(meta_path))
               jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
             else
               NULL

  outer_folds  <- as.integer(meta$outer_folds)
  inner_folds  <- as.integer(meta$inner_folds)
  feature_cols <- setdiff(meta$features, "row_id")
  if (!is.null(drop_cols_after_impute))
    feature_cols <- setdiff(feature_cols, drop_cols_after_impute)

  time_var  <- meta$time_col
  event_var <- meta$event_col

  meta     <- get_cr_model(model_key)$info()
  has_grid <- isTRUE(meta$needs_tuning)

  if (has_grid) {
    if (is.null(grid))
      stop("Model '", model_key, "' requires tuning; please supply a `grid`.")
    grid_df   <- tidyr::expand_grid(!!!grid)
    grid_list <- lapply(seq_len(nrow(grid_df)),
                        function(i) as.list(grid_df[i, , drop = FALSE]))
    if (verbose) {
      cat("Model:", meta$name, "\n")
      cat("Outer folds:", outer_folds, "| Inner folds:", inner_folds, "\n")
      cat("Has tuning: TRUE | #grid cfgs:", length(grid_list), "\n")
    }
  } else {
    if (verbose) {
      cat("Model:", meta$name, "\n")
      cat("Outer folds:", outer_folds, "\n")
    }
  }

  results <- vector("list", outer_folds)

  for (v in seq_len(outer_folds)) {
    fold_dir   <- file.path(out_dir, sprintf("outer_%d", v))
    outer_data <- .load_outer_data(fold_dir, cr_meta)
    m_outer    <- length(outer_data)

    train  <- outer_data[[1L]]$train
    test   <- outer_data[[1L]]$test
    causes <- sort(unique(train$event[train$event != 0L]))

    set.seed(seed + v)

    # ---- no tuning ----
    if (!has_grid) {
      cif_list <- vector("list", m_outer)
      for (m in seq_len(m_outer)) {
        tr_m <- outer_data[[m]]$train
        te_m <- outer_data[[m]]$test
        fit  <- fit_cr_model(model_key,
                               obj  = cr_data(tr_m, time_var = time_var,
                                              event_var = event_var),
                               args = list())
        cif_list[[m]] <- predict_cif(fit,
          newdata   = cr_data(te_m, time_var = time_var, event_var = event_var),
          time_grid = times)
      }
      pred    <- pool_cifs_mean(cif_list)
      cr_test <- cr_data(test, time_var = time_var, event_var = event_var)
      perf    <- compute_metrics(cr_test,
                                 cif           = list(cif = pred, time_grid = times),
                                 pred_horizons = times,
                                 metrics       = c("Brier", "IBS"),
                                 collapse_as_df = FALSE)

      store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
      save_cif_r(store_dir, cif = pred, times = times,
                 row_ids = rownames(test), causes = causes)

      results[[v]] <- list(
        fold      = v,
        model_key = model_key,
        times     = times,
        metrics   = list(causes    = causes,
                         cause_bs  = perf$Brier,
                         cause_ibs = perf$IBS)
      )
      if (verbose) cat("[outer", v, "] done (no tuning) \n")

    # ---- with tuning ----
    } else {
      inner_loaded <- .load_inner_data(fold_dir, cr_meta,
                                      outer_train_df = train)
      inner_dirs   <- inner_loaded$inner_dirs
      inner_data   <- inner_loaded$inner_data

      cfg_scores <- sapply(seq_along(grid_list), function(gi) {
        cfg <- grid_list[[gi]]
        if (verbose)
          cat("[outer", v, "] grid", gi, "/", length(grid_list), ":",
              paste(names(cfg), cfg, collapse = " "), "\n")

        inner_vals <- matrix(NA_real_,
                             nrow = length(causes), ncol = inner_folds)

        for (j in seq_along(inner_dirs)) {
          jj       <- inner_data[[j]]
          m_inner  <- length(jj)
          cif_list <- vector("list", m_inner)

          for (m in seq_len(m_inner)) {
            tr_m   <- jj[[m]]$train
            va_m   <- jj[[m]]$val
            fit_in <- fit_cr_model(model_key,
                                   obj  = cr_data(tr_m, time_var = time_var,
                                                  event_var = event_var),
                                   args = cfg)
            cif_list[[m]] <- predict_cif(fit_in,
              newdata   = cr_data(va_m, time_var = time_var, event_var = event_var),
              time_grid = times)
          }

          cif_in   <- pool_cifs_mean(cif_list)
          val_ref  <- jj[[1L]]$val
          cr_val   <- cr_data(val_ref, time_var = time_var, event_var = event_var)
          perf_in  <- compute_metrics(cr_val,
                                      cif           = list(cif = cif_in, time_grid = times),
                                      pred_horizons = times,
                                      metrics       = "IBS",
                                      collapse_as_df = FALSE)

          for (ci in seq_along(causes)) {
            kk <- causes[ci]
            inner_vals[ci, j] <- utils::tail(
              perf_in$IBS[[paste0("cause_", kk)]], 1L
            )
          }
        }
        inner_vals
      })

      ncause <- length(causes)
      ngrids <- length(grid_list)

      cfg_mat <- array(as.vector(cfg_scores),
                       dim = c(ncause, inner_folds, ngrids),
                       dimnames = list(
                         cause = paste0("cause_", causes),
                         fold  = paste0("fold_",  seq_len(inner_folds)),
                         grid  = paste0("cfg_",   seq_len(ngrids))
                       ))
      cfg_av  <- apply(cfg_mat, c(1L, 3L), mean, na.rm = TRUE)
      best_idx <- apply(cfg_av, 1L, which.min)
      best_cfg <- mapply(function(i) grid_list[[i]], best_idx,
                         SIMPLIFY = FALSE)

      pred <- array(NA_real_, dim = c(nrow(test), ncause, length(times)))

      frgp_details <- identical(model_key, "FGRP")
      if (frgp_details) {
        store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
        dir.create(store_dir, recursive = TRUE, showWarnings = FALSE)
      }

      for (ci in seq_along(causes)) {
        kk    <- causes[ci]
        cfg_i <- best_cfg[[paste0("cause_", kk)]]
        cif_list  <- vector("list", m_outer)
        coef_list <- if (frgp_details) vector("list", m_outer) else NULL
        fit_list  <- if (frgp_details) vector("list", m_outer) else NULL

        for (mii in seq_len(m_outer)) {
          tr_k  <- outer_data[[mii]]$train
          te_k  <- outer_data[[mii]]$test
          fit_i <- fit_cr_model(model_key,
                                 obj  = cr_data(tr_k, time_var = time_var,
                                                event_var = event_var),
                                 args = cfg_i)
          cif_list[[mii]] <- predict_cif(fit_i,
            newdata   = cr_data(te_k, time_var = time_var, event_var = event_var),
            time_grid = times)
          if (frgp_details) {
            fit_list[[mii]]  <- fit_i
            coef_list[[mii]] <- extract_fgrp_coefs(fit_i, tol = 1e-8)
          }
        }

        pred_i          <- pool_cifs_mean(cif_list)
        pred[, ci, ]    <- pred_i[, ci, ]

        if (verbose)
          cat("[outer", v, "] best cfg idx:",
              best_idx[paste0("cause_", kk)],
              " avg.score (across inner folds):",
              sprintf("%.6f", cfg_av[paste0("cause_", kk),
                                     paste0("cfg_", best_idx[paste0("cause_", kk)])]),
              " cause:", kk, "\n")

        if (frgp_details) {
          store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
          saveRDS(cfg_i,
                  file.path(store_dir, sprintf("best_cfg_cause_%s.rds", kk)),
                  compress = "xz")
          saveRDS(fit_list,
                  file.path(store_dir, sprintf("best_fit_cause_%s.rds", kk)),
                  compress = "xz")
          coef_df <- do.call(rbind, lapply(seq_along(coef_list), function(i)
            cbind(imputation = i, coef_list[[i]])))
          utils::write.csv(coef_df,
                           file.path(store_dir,
                                     sprintf("best_coef_cause_%s.csv", kk)),
                           row.names = FALSE)
        }
      }

      store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
      save_cif_r(store_dir, cif = pred, times = times,
                 row_ids = rownames(test), causes = causes)

      cr_test <- cr_data(test, time_var = time_var, event_var = event_var)
      perf    <- compute_metrics(cr_test,
                                 cif           = list(cif = pred, time_grid = times),
                                 pred_horizons = times,
                                 metrics       = "IBS",
                                 collapse_as_df = FALSE)

      results[[v]] <- list(
        fold      = v,
        model_key = model_key,
        times     = times,
        best_cfg  = best_cfg,
        metrics   = list(causes    = causes,
                         cause_ibs = perf$IBS)
      )
    }
  }
  results
}


