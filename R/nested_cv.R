#' @title Nested Cross-Validation for Competing Risks Models
#' @description The core engine for running nested cross-validation over a
#'   pre-built benchmark directory created by [create_nested_splits()].
#' @name nested_cv
NULL

# ---- internal helpers --------------------------------------------------------

#' @noRd
read_bench_parquet <- function(path, schema) {
  df <- apply_data_types(
    as.data.frame(arrow::read_parquet(path)),
    schema = schema
  )
  if (!("row_id" %in% names(df)))
    stop("row_id column not found in: ", path, call. = FALSE)
  rownames(df) <- as.character(df$row_id)
  df$row_id    <- NULL
  df
}

#' @noRd
detect_m_from_dir <- function(dir) {
  f  <- list.files(dir, pattern = "^train_imp_[0-9]+\\.parquet$",
                   full.names = FALSE)
  if (!length(f)) return(1L)
  ks <- sort(unique(as.integer(
    sub("^train_imp_([0-9]+)\\.parquet$", "\\1", f)
  )))
  ks <- ks[!is.na(ks)]
  if (!length(ks) || ks[1L] != 1L)
    stop("Invalid train_imp_*.parquet numbering in: ", dir, call. = FALSE)
  length(ks)
}

#' @noRd
load_outer_data <- function(fold_dir, schema) {
  if (file.exists(file.path(fold_dir, "train_imp_1.parquet")) &&
      file.exists(file.path(fold_dir, "test_imp_1.parquet"))) {
    mm  <- detect_m_from_dir(fold_dir)
    out <- vector("list", mm)
    for (k in seq_len(mm)) {
      trp <- file.path(fold_dir, sprintf("train_imp_%d.parquet", k))
      tep <- file.path(fold_dir, sprintf("test_imp_%d.parquet",  k))
      out[[k]] <- list(train = read_bench_parquet(trp, schema),
                       test  = read_bench_parquet(tep, schema))
    }
    return(out)
  }
  tr_path <- file.path(fold_dir, "train.parquet")
  te_path <- file.path(fold_dir, "test.parquet")
  list(list(train = read_bench_parquet(tr_path, schema),
            test  = read_bench_parquet(te_path, schema)))
}

#' @noRd
load_inner_data <- function(fold_dir, schema, outer_train_df) {
  inner_dir  <- file.path(fold_dir, "inner")
  inner_dirs <- if (dir.exists(inner_dir)) {
    dd <- list.dirs(inner_dir, recursive = FALSE, full.names = TRUE)
    dd[grepl("fold_[0-9]+$", basename(dd))]
  } else character(0L)

  if (!length(inner_dirs))
    return(list(inner_dirs = inner_dirs, inner_data = list()))

  inner_data <- vector("list", length(inner_dirs))

  for (j in seq_along(inner_dirs)) {
    idir <- inner_dirs[j]
    if (file.exists(file.path(idir, "train_imp_1.parquet")) &&
        file.exists(file.path(idir, "val_imp_1.parquet"))) {
      mm <- detect_m_from_dir(idir)
      jj <- vector("list", mm)
      for (k in seq_len(mm)) {
        trp <- file.path(idir, sprintf("train_imp_%d.parquet", k))
        vap <- file.path(idir, sprintf("val_imp_%d.parquet",   k))
        jj[[k]] <- list(train = read_bench_parquet(trp, schema),
                        val   = read_bench_parquet(vap, schema))
      }
      inner_data[[j]] <- jj
      next
    }
    idx_path <- file.path(idir, "indices.json")
    idx      <- jsonlite::fromJSON(idx_path, simplifyVector = TRUE)
    tr_ids   <- as.character(idx$train_row_id)
    va_ids   <- as.character(idx$val_row_id)
    inner_data[[j]] <- list(list(
      train = outer_train_df[tr_ids, , drop = FALSE],
      val   = outer_train_df[va_ids, , drop = FALSE]
    ))
  }

  list(inner_dirs = inner_dirs, inner_data = inner_data)
}


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
  schema_path   <- file.path(out_dir, "data_types.json")

  man    <- jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
  times  <- sort(unique(as.numeric(
    jsonlite::fromJSON(times_path, simplifyVector = TRUE)
  )))
  schema <- jsonlite::fromJSON(schema_path, simplifyVector = TRUE)

  outer_folds  <- as.integer(man$outer_folds)
  inner_folds  <- as.integer(man$inner_folds)
  feature_cols <- setdiff(man$features, "row_id")
  if (!is.null(drop_cols_after_impute))
    feature_cols <- setdiff(feature_cols, drop_cols_after_impute)

  # Time and event variable names — standardised by prepare_data()
  time_var  <- "time"
  event_var <- "event"

  mdl      <- get_cr_model(model_key)
  meta     <- mdl$info()
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
    outer_data <- load_outer_data(fold_dir, schema)
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
        fit  <- mdl$fit(obj  = cr_data(tr_m, time_var = time_var,
                                         event_var = event_var),
                        args = list())
        cif_list[[m]] <- mdl$predict_cif(fit,
          newdata   = cr_data(te_m, time_var = time_var, event_var = event_var),
          time_grid = times)
      }
      pred  <- pool_cifs_mean(cif_list)
      perf  <- score_from_cifs(out = pred, test = test,
                               times = times, causes = causes)

      store_dir <- build_store_paths_r(out_dir, model = model_key, fold = v)
      save_cif_r(store_dir, cif = pred, times = times,
                 row_ids = rownames(test), causes = causes)

      results[[v]] <- list(
        fold      = v,
        model_key = model_key,
        times     = times,
        metrics   = list(causes   = fit$causes,
                         cause_bs  = perf$cause_bs,
                         cause_ibs = perf$cause_ibs)
      )
      if (verbose) cat("[outer", v, "] done (no tuning) \n")

    # ---- with tuning ----
    } else {
      inner_loaded <- load_inner_data(fold_dir, schema,
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
            fit_in <- mdl$fit(obj  = cr_data(tr_m, time_var = time_var,
                                               event_var = event_var),
                              args = cfg)
            cif_list[[m]] <- mdl$predict_cif(fit_in,
              newdata   = cr_data(va_m, time_var = time_var, event_var = event_var),
              time_grid = times)
          }

          cif_in  <- pool_cifs_mean(cif_list)
          val_ref <- jj[[1L]]$val
          perf_in <- score_from_cifs(out = cif_in, test = val_ref,
                                     times = times, causes = causes)

          for (ci in seq_along(causes)) {
            kk <- causes[ci]
            inner_vals[ci, j] <- utils::tail(
              perf_in$cause_ibs[[kk]], 1L
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
          fit_i <- mdl$fit(obj  = cr_data(tr_k, time_var = time_var,
                                            event_var = event_var),
                           args = cfg_i)
          cif_list[[mii]] <- mdl$predict_cif(fit_i,
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

      perf <- score_from_cifs(out = pred, test = test,
                              times = times, causes = causes)

      results[[v]] <- list(
        fold      = v,
        model_key = model_key,
        times     = times,
        best_cfg  = best_cfg,
        metrics   = list(causes    = causes,
                         cause_ibs = perf$cause_ibs)
      )
    }
  }
  results
}


