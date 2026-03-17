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
#' @param cif_time_grid Numeric vector of time points at which CIF predictions
#'   are evaluated.
#' @param grid      For tunable models: a named list of hyperparameter
#'   sequences to be expanded via `tidyr::expand_grid()`.  Ignored for
#'   models with `needs_tuning = FALSE`.
#' @param seed      Base random seed (outer fold `v` uses `seed + v`).
#' @param verbose   Print progress messages.
#' @param drop_cols_after_impute Optional character vector of feature columns
#'   to exclude from the model even after they have been used in imputation.
#'
#' @return A list of length `outer_folds`.  Each element is a named list
#'   containing `fold`, `model_key`, `cif_time_grid`, `best_cfg` (tunable models
#'   only), and `metrics` (`causes`, `cause_bs`, `cause_ibs`).
#' @export
nested_cv_from_bench <- function(out_dir                = "../BenchResults",
                                 model_key,
                                 cif_time_grid,
                                 grid                   = NULL,
                                 seed                   = 123L,
                                 verbose                = TRUE,
                                 drop_cols_after_impute = NULL) {
  
  # Load CV metadata
  cv_meta_path <- file.path(out_dir, "cv_splits_metadata.json")
  cv_meta <- jsonlite::fromJSON(cv_meta_path, simplifyVector = TRUE)
  
  # Extract CV metadata - consider removing drop_cols_after_impute (this would be handled by the imputation function)
  outer_folds  <- as.integer(cv_meta$outer_folds)
  inner_folds  <- as.integer(cv_meta$inner_folds)
  feature_cols <- setdiff(cv_meta$var_names,
                          c(cv_meta$time_var, cv_meta$event_var, cv_meta$id_var))
  if (!is.null(drop_cols_after_impute))
    feature_cols <- setdiff(feature_cols, drop_cols_after_impute)
  
  time_var  <- cv_meta$time_var
  event_var <- cv_meta$event_var
  id_var    <- cv_meta$id_var
  causes    <- as.integer(cv_meta$causes)
  
  # Get information about the model and whether it requires tuning
  model_info <- get_cr_model(model_key)$info()
  has_grid   <- isTRUE(model_info$needs_tuning)
  
  # Define grid for tuning parameters and show progress messages
  if (has_grid) {
    if (is.null(grid))
      stop("Model '", model_key, "' requires tuning; please supply a `grid`.")
    grid_df   <- tidyr::expand_grid(!!!grid)
    grid_list <- lapply(seq_len(nrow(grid_df)),
                        function(i) as.list(grid_df[i, , drop = FALSE]))
    if (verbose) {
      cat("Model:", model_info$name, "\n")
      cat("Outer folds:", outer_folds, "| Inner folds:", inner_folds, "\n")
      cat("Has tuning: TRUE | #grid cfgs:", length(grid_list), "\n")
    }
  } else {
    if (verbose) {
      cat("Model:", model_info$name, "\n")
      cat("Outer folds:", outer_folds, "\n")
    }
  }
  
  # Create output object to store results
  results <- vector("list", outer_folds)
  
  # Loop over the outer folds
  for (v in seq_len(outer_folds)) {
    
    # Read the outer fold data (train/test splits, possibly with multiple imputations)
    fold_dir   <- file.path(out_dir, sprintf("outer_%d", v))
    outer_data <- .load_outer_data(fold_dir, cv_meta)
    # Note: m_outer > 1 if multiple imputations were stored for the outer fold; otherwise m_outer = 1
    m_outer    <- length(outer_data)
    
    train  <- outer_data[[1]]$train
    test   <- outer_data[[1]]$test
    
    set.seed(seed + v)
    
    # ---- no tuning ----
    if (!has_grid) {
      cif_list <- vector("list", m_outer)
      for (m in seq_len(m_outer)) {
        tr_m <- outer_data[[m]]$train
        te_m <- outer_data[[m]]$test
        # Fit a model using the training set
        fit  <- fit_cr_model(model_key,
                             obj  = cr_data(tr_m, time_var = time_var,
                                            event_var = event_var,
                                            id_var = id_var),
                             args = list())
        # Predict the CIF on the test set at the specified time points and store in list
        cif_list[[m]] <- predict_cif(
          fit,
          newdata   = cr_data(te_m, time_var = time_var, event_var = event_var, id_var = id_var),
          time_grid = cif_time_grid)
      }
      # Pool CIF predictions across imputations (no-op when m_outer = 1);
      # pred is a predict_cif-style list with $cif, $time_grid, $model_key, $ids
      pred <- .pool_cifs_mean(cif_list)
      
      store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
      save_cif_r(store_dir, cif = pred$cif, cif_time_grid = pred$time_grid,
                 row_ids = rownames(test), causes = causes)
      
      results[[v]] <- list(
        fold          = v,
        model_key     = model_key,
        cif_time_grid = cif_time_grid
      )
      if (verbose) cat("[outer", v, "] done (no tuning) \n")
      
      # ---- with tuning ----
    } else {
      # Read the inner fold data (train/val splits within the outer training set)
      inner_loaded <- .load_inner_data(fold_dir, cv_meta,
                                       outer_train_df = train)
      inner_dirs   <- inner_loaded$inner_dirs
      inner_data   <- inner_loaded$inner_data
      
      # Evaluate each hyperparameter configuration across all inner folds.
      # cfg_scores is a matrix [ncause x inner_folds] per grid config, returned
      # as a list-column by sapply — reshaped into a 3D array below.
      cfg_scores <- sapply(seq_along(grid_list), function(gi) {
        cfg <- grid_list[[gi]]
        if (verbose)
          cat("[outer", v, "] grid", gi, "/", length(grid_list), ":",
              paste(names(cfg), cfg, collapse = " "), "\n")
        
        # Accumulate IBS per cause per inner fold for this config
        inner_vals <- matrix(NA_real_,
                             nrow = length(causes), ncol = inner_folds)
        
        for (j in seq_along(inner_dirs)) {
          jj       <- inner_data[[j]]
          # Note: m_inner > 1 if multiple imputations are present; otherwise m_inner = 1
          m_inner  <- length(jj)
          cif_list <- vector("list", m_inner)
          
          for (m in seq_len(m_inner)) {
            tr_m   <- jj[[m]]$train
            va_m   <- jj[[m]]$val
            # Fit model on inner training set with candidate hyperparameters
            fit_in <- fit_cr_model(model_key,
                                   obj  = cr_data(tr_m, time_var = time_var,
                                                  event_var = event_var,
                                                  id_var = id_var),
                                   args = cfg)
            # Predict CIF on inner validation set and store in list
            cif_list[[m]] <- predict_cif(fit_in,
                                         newdata   = cr_data(va_m, time_var = time_var, event_var = event_var, id_var = id_var),
                                         time_grid = cif_time_grid)
          }
          
          # Pool CIF predictions across imputations (no-op when m_inner = 1)
          cif_in  <- .pool_cifs_mean(cif_list)
          # Use true (unimputed) outcomes from the first imputation's val set
          val_ref <- jj[[1L]]$val
          cr_val  <- cr_data(val_ref, time_var = time_var, event_var = event_var, id_var = id_var)
          # Compute IBS on the pooled CIF for this inner fold
          perf_in <- compute_metrics(cr_val,
                                     cif           = cif_in,
                                     pred_horizons = cif_time_grid,
                                     metrics       = "IBS",
                                     collapse_as_df = FALSE)
          
          # Store the cumulative IBS (last value = IBS up to max time) per cause
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
      
      # Reshape cfg_scores into a [ncause x inner_folds x ngrids] array and
      # average across inner folds to get mean IBS per cause per config
      cfg_mat <- array(as.vector(cfg_scores),
                       dim = c(ncause, inner_folds, ngrids),
                       dimnames = list(
                         cause = paste0("cause_", causes),
                         fold  = paste0("fold_",  seq_len(inner_folds)),
                         grid  = paste0("cfg_",   seq_len(ngrids))
                       ))
      cfg_av   <- apply(cfg_mat, c(1L, 3L), mean, na.rm = TRUE)
      # Select the config with the lowest average IBS per cause independently
      best_idx <- apply(cfg_av, 1L, which.min)
      best_cfg <- mapply(function(i) grid_list[[i]], best_idx,
                         SIMPLIFY = FALSE)
      
      # Initialise the raw prediction array for the outer test set;
      # causes are filled in separately below and pooled into a single array
      pred <- array(NA_real_, dim = c(nrow(test), ncause, length(cif_time_grid)))
      
      frgp_details <- identical(model_key, "FGRP")
      if (frgp_details) {
        store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
        dir.create(store_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Refit using the best config per cause on the full outer training set
      for (ci in seq_along(causes)) {
        kk        <- causes[ci]
        cfg_i     <- best_cfg[[paste0("cause_", kk)]]
        cif_list  <- vector("list", m_outer)
        coef_list <- if (frgp_details) vector("list", m_outer) else NULL
        fit_list  <- if (frgp_details) vector("list", m_outer) else NULL
        
        for (mii in seq_len(m_outer)) {
          tr_k  <- outer_data[[mii]]$train
          te_k  <- outer_data[[mii]]$test
          # Fit model on full outer training set using the best config for this cause
          fit_i <- fit_cr_model(model_key,
                                obj  = cr_data(tr_k, time_var = time_var,
                                               event_var = event_var,
                                               id_var = id_var),
                                args = cfg_i)
          # Predict CIF on outer test set and store in list
          cif_list[[mii]] <- predict_cif(fit_i,
                                         newdata   = cr_data(te_k, time_var = time_var, event_var = event_var, id_var = id_var),
                                         time_grid = cif_time_grid)
          if (frgp_details) {
            fit_list[[mii]]  <- fit_i
            coef_list[[mii]] <- extract_fgrp_coefs(fit_i, tol = 1e-8)
          }
        }
        
        # Pool CIF predictions across imputations and slot into the combined array
        pred_i       <- .pool_cifs_mean(cif_list)
        pred[, ci, ] <- pred_i$cif[, ci, ]
        
        if (verbose)
          cat("[outer", v, "] best cfg idx:",
              best_idx[paste0("cause_", kk)],
              " avg.score (across inner folds):",
              sprintf("%.6f", cfg_av[paste0("cause_", kk),
                                     paste0("cfg_", best_idx[paste0("cause_", kk)])]),
              " cause:", kk, "\n")
        
        # FGRP-specific: save best config, fitted models, and coefficients to disk
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
      save_cif_r(store_dir, cif = pred, cif_time_grid = cif_time_grid,
                 row_ids = rownames(test), causes = causes)
      
      results[[v]] <- list(
        fold          = v,
        model_key     = model_key,
        cif_time_grid = cif_time_grid,
        best_cfg      = best_cfg
      )
    }
  }
  results
}