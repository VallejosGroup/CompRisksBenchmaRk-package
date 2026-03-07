#' @title Data Utilities for Competing Risks Benchmarking
#' @description Functions for preparing data, creating nested cross-validation
#'   splits, imputing missing values, and miscellaneous helpers.
#' @name utils
NULL

#' Prepare a data frame for benchmarking
#'
#' Renames the time and event columns to `"time"` and `"event"`, adds a small
#' offset when the minimum observed time is zero (to avoid metric failures),
#' and sorts rows by time.
#'
#' @param df A data frame.
#' @param time_col Name of the column containing observed times.
#' @param event_col Name of the column containing event/status codes
#'   (0 = censored).
#'
#' @return The modified data frame sorted by time.
#' @export
prepare_data <- function(df, time_col, event_col) {
  nm <- names(df)
  nm[nm == time_col]  <- "time"
  nm[nm == event_col] <- "event"
  names(df) <- nm
  
  min_time <- suppressWarnings(min(df$time, na.rm = TRUE))
  if (is.finite(min_time) && min_time == 0) {
    df$time <- df$time + 0.01
  }
  
  df[order(df$time), , drop = FALSE]
}


#' Infer column data types from a data frame
#'
#' Returns a schema list compatible with [apply_data_types()].
#'
#' @param df A data frame.
#' @param feature_cols Character vector of feature column names.
#' @param core_cols Character vector of mandatory non-feature columns
#'   (default `c("time","event","row_id")`).
#'
#' @return A named list with elements `version`, `core_cols`,
#'   `feature_cols`, and `types`.
#' @export
infer_data_types <- function(df, feature_cols,
                             core_cols = c("time", "event", "row_id")) {
  types <- stats::setNames(
    vector("list", length(core_cols) + length(feature_cols)),
    c(core_cols, feature_cols)
  )
  
  types[["time"]]   <- "numeric"
  types[["event"]]  <- "integer"
  types[["row_id"]] <- "character"
  
  for (cn in feature_cols) {
    if (!cn %in% names(df)) next
    x <- df[[cn]]
    types[[cn]] <- if (is.logical(x))   "logical"
    else if (is.factor(x))    "factor"
    else if (is.character(x)) "character"
    else if (is.integer(x))   "integer"
    else if (is.numeric(x))   "numeric"
    else                      "character"
  }
  
  list(version = 1L, core_cols = core_cols,
       feature_cols = feature_cols, types = as.list(types))
}


#' Create nested cross-validation splits and write them to disk
#'
#' Produces a directory hierarchy suitable for [nested_cv_from_bench()].
#' Outer and inner folds are stratified by event status using
#' [rsample::vfold_cv()].  Fold data are stored as Parquet files; inner
#' fold indices are stored as JSON.
#'
#' @param df A data frame containing at minimum time, event, and covariate
#'   columns.
#' @param time_col Name of the time column.
#' @param event_col Name of the event/status column.
#' @param id_col Optional name of an existing ID column.  If `NULL`, row
#'   numbers are used.
#' @param out_dir Path to the output directory (created if needed).
#' @param outer_folds Number of outer folds (default 5).
#' @param inner_folds Number of inner folds (default 3).
#' @param seed Random seed for reproducibility.
#' @param times Numeric vector of evaluation time points written to
#'   `times.json`.
#'
#' @return Invisibly `TRUE`.
#' @export
create_nested_splits <- function(df,
                                 time_col,
                                 event_col,
                                 id_col       = NULL,
                                 out_dir      = NULL,
                                 outer_folds  = 5,
                                 inner_folds  = 3,
                                 seed         = 123,
                                 times) {
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  df <- prepare_data(df = df, time_col = time_col, event_col = event_col)
  
  if (!is.null(id_col)) {
    df$row_id <- as.character(df[[id_col]])
    df[[id_col]] <- NULL
  } else {
    df$row_id <- as.character(seq_len(nrow(df)))
  }
  
  core_cols    <- c("time", "event", "row_id")
  feature_cols <- setdiff(names(df), core_cols)
  
  schema <- infer_data_types(df, feature_cols = feature_cols,
                             core_cols = core_cols)
  jsonlite::write_json(schema,
                       file.path(out_dir, "data_types.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  
  for_python <- list(
    version     = 1,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    features    = feature_cols,
    time_col    = "time",
    event_col   = "event",
    row_id      = "row_id"
  )
  jsonlite::write_json(for_python,
                       file.path(out_dir, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  jsonlite::write_json(times,
                       file.path(out_dir, "times.json"),
                       auto_unbox = TRUE)
  
  set.seed(seed)
  outer <- rsample::vfold_cv(df, v = outer_folds, strata = "event")
  
  for (v in seq_len(outer_folds)) {
    split_v <- outer$splits[[v]]
    train   <- rsample::analysis(split_v)
    test    <- rsample::assessment(split_v)
    
    out_v_dir <- file.path(out_dir, sprintf("outer_%d", v))
    dir.create(out_v_dir, recursive = TRUE, showWarnings = FALSE)
    
    arrow::write_parquet(
      dplyr::select(train, dplyr::all_of(c("time", "event", feature_cols, "row_id"))),
      file.path(out_v_dir, "train.parquet")
    )
    arrow::write_parquet(
      dplyr::select(test,  dplyr::all_of(c("time", "event", feature_cols, "row_id"))),
      file.path(out_v_dir, "test.parquet")
    )
    
    set.seed(seed + v)
    inner     <- rsample::vfold_cv(train, v = inner_folds, strata = "event")
    inner_dir <- file.path(out_v_dir, "inner")
    dir.create(inner_dir, showWarnings = FALSE)
    
    for (j in seq_len(inner_folds)) {
      split_j    <- inner$splits[[j]]
      train_in   <- rsample::analysis(split_j)
      val_in     <- rsample::assessment(split_j)
      fold_j_dir <- file.path(inner_dir, sprintf("fold_%d", j))
      dir.create(fold_j_dir, showWarnings = FALSE)
      
      idx_payload <- list(
        train_row_id = train_in[["row_id"]],
        val_row_id   = val_in[["row_id"]]
      )
      jsonlite::write_json(idx_payload,
                           file.path(fold_j_dir, "indices.json"),
                           auto_unbox = TRUE)
    }
  }
  
  invisible(TRUE)
}


#' Apply stored data types to a data frame
#'
#' Coerces columns of `df` to the types recorded in `schema` (as produced
#' by [infer_data_types()]).
#'
#' @param df A data frame.
#' @param schema A schema list as returned by [infer_data_types()].
#'
#' @return The coerced data frame.
#' @export
apply_data_types <- function(df, schema) {
  out   <- df
  types <- schema$types
  
  if ("time"   %in% names(out)) out$time   <- as.numeric(out$time)
  if ("event"  %in% names(out)) out$event  <- as.integer(out$event)
  if ("row_id" %in% names(out)) out$row_id <- as.character(out$row_id)
  
  for (cn in names(types)) {
    if (!cn %in% names(out)) next
    if (cn %in% c("time", "event", "row_id")) next
    typ <- types[[cn]]
    out[[cn]] <- switch(typ,
                        numeric   = as.numeric(out[[cn]]),
                        integer   = as.integer(out[[cn]]),
                        logical   = as.logical(out[[cn]]),
                        factor    = as.factor(out[[cn]]),
                        as.character(out[[cn]])
    )
  }
  out
}


#' Detect missing values in train/test splits
#'
#' @param train Training data frame.
#' @param test  Test data frame.
#' @param cols  Character vector of column names to check.
#'
#' @return `TRUE` if any NA is found, `FALSE` otherwise.
#' @keywords internal
split_detect_NA <- function(train, test, cols) {
  anyNA(train[, cols, drop = FALSE]) ||
    anyNA(test[, cols, drop = FALSE])
}


# ---- Imputation helpers ----

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA)
  tab <- table(x)
  names(tab)[which.max(tab)]
}

impute_median_mode_single <- function(train, test, feature_cols,
                                      drop_cols_after_impute) {
  tr <- train
  ts <- test
  message("Imputation: median/mode")
  for (feat in feature_cols) {
    f <- tr[[feat]]
    if (is.numeric(f) || is.integer(f)) {
      med <- stats::median(f, na.rm = TRUE)
      tr[[feat]][is.na(tr[[feat]])] <- med
      ts[[feat]][is.na(ts[[feat]])] <- med
    } else if (is.factor(f) || is.logical(f)) {
      fill <- mode_value(f)
      tr[[feat]][is.na(tr[[feat]])] <- fill
      ts[[feat]][is.na(ts[[feat]])] <- fill
    }
  }
  if (!is.null(drop_cols_after_impute)) {
    tr <- tr[, !(names(tr) %in% drop_cols_after_impute)]
    ts <- ts[, !(names(ts) %in% drop_cols_after_impute)]
  }
  list(train = tr, test = ts)
}


stop_if_any_na <- function(x, name = deparse(substitute(x))) {
  if (anyNA(x)) {
    na_cells <- sum(is.na(x))
    na_rows  <- sum(rowSums(is.na(x)) > 0)
    na_cols  <- sum(colSums(is.na(x)) > 0)
    top_cols <- sort(colSums(is.na(x)), decreasing = TRUE)
    top_cols <- top_cols[top_cols > 0]
    stop(
      sprintf(
        "%s contains NA after imputation: %d NA cells across %d rows and %d cols. Top NA cols: %s",
        name, na_cells, na_rows, na_cols,
        paste(sprintf("%s=%d", names(top_cols), as.integer(top_cols)),
              collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


impute_mice <- function(train, test, feature_cols, m, maxit, seed,
                        drop_cols_after_impute) {
  train_x <- train[, feature_cols, drop = FALSE]
  test_x  <- test[,  feature_cols, drop = FALSE]
  
  meth  <- mice::make.method(train_x)
  predM <- mice::make.predictorMatrix(train_x)
  
  imp_train <- mice::mice(train_x, m = m, method = meth,
                          predictorMatrix = predM, maxit = maxit,
                          seed = seed, printFlag = FALSE, ridge = 1e-6)
  imp_test  <- mice::mice.mids(imp_train, newdata = test_x,
                               printFlag = FALSE)
  
  train_x_imp <- mice::complete(imp_train, action = "all")
  test_x_imp  <- mice::complete(imp_test,  action = "all")
  
  keep_cols    <- setdiff(names(train), feature_cols)
  imputations  <- vector("list", m)
  
  for (k in seq_len(m)) {
    tr_k <- cbind(train[, keep_cols, drop = FALSE], train_x_imp[[k]])
    te_k <- cbind(test[,  keep_cols, drop = FALSE], test_x_imp[[k]])
    
    if (!is.null(drop_cols_after_impute)) {
      tr_k <- tr_k[, !(names(tr_k) %in% drop_cols_after_impute), drop = FALSE]
      te_k <- te_k[, !(names(te_k) %in% drop_cols_after_impute), drop = FALSE]
    }
    stop_if_any_na(tr_k, sprintf("MICE train (k=%d)", k))
    stop_if_any_na(te_k, sprintf("MICE test  (k=%d)", k))
    imputations[[k]] <- list(train = tr_k, test = te_k)
  }
  list(m = m, imputations = imputations)
}


#' Impute missing values in a train/test split
#'
#' Handles three strategies: `"none"` (passthrough), `"median_mode"`
#' (numeric columns filled with training median, factor/logical columns
#' filled with training mode), and `"mice"` (multiple imputation via the
#' **mice** package).
#'
#' @param train Training data frame.
#' @param test  Test (or validation) data frame.
#' @param feature_cols Character vector of feature column names to impute.
#' @param method One of `"none"`, `"median_mode"`, or `"mice"`.
#' @param m Number of imputed datasets (only used when `method = "mice"`).
#' @param maxit Number of MICE iterations (only used when `method = "mice"`).
#' @param seed Integer seed for MICE.
#' @param drop_cols_after_impute Optional character vector of column names
#'   to drop after imputation.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{m}{Number of imputed datasets.}
#'     \item{imputations}{A list of length `m`, each element being a list
#'       with `$train` and `$test` data frames.}
#'   }
#' @export
impute_split <- function(train, test, feature_cols,
                         method = c("none", "median_mode", "mice"),
                         m      = 1L,
                         maxit  = 100L,
                         seed   = 1L,
                         drop_cols_after_impute = NULL) {
  method <- match.arg(method)
  
  if (method == "mice") {
    cat("Default number of datasets m: ", m, ", and iterations: ", maxit,
        " seed: ", seed, "\n")
  }
  
  if (!split_detect_NA(train, test, feature_cols) || method == "none") {
    return(list(m = 1L, imputations = list(list(train = train, test = test))))
  }
  
  if (method == "median_mode") {
    imp <- impute_median_mode_single(
      train = train, test = test,
      feature_cols = feature_cols,
      drop_cols_after_impute = drop_cols_after_impute
    )
    return(list(m = 1L,
                imputations = list(list(train = imp$train, test = imp$test))))
  }
  
  if (method == "mice") {
    return(impute_mice(
      train = train, test = test,
      feature_cols = feature_cols, m = m, maxit = maxit, seed = seed,
      drop_cols_after_impute = drop_cols_after_impute
    ))
  }
  
  stop("Unsupported method: ", method)
}


#' Pool a list of CIF arrays by taking the element-wise mean
#'
#' Used to combine predictions across multiple imputed datasets.
#'
#' @param cif_list A list of 3-D numeric arrays, each of dimension
#'   `[n, K, Tm]`.
#'
#' @return A single 3-D array of the same dimensions.
#' @export
pool_cifs_mean <- function(cif_list) {
  m      <- length(cif_list)
  pooled <- cif_list[[1]]
  if (m > 1L)
    for (i in 2L:m) pooled <- pooled + cif_list[[i]]
  pooled / m
}


#' Build a numeric design matrix aligned between train and new data
#'
#' Expands factor columns, drops zero-variance columns, enforces full
#' column rank, and ensures `x_new` has exactly the same columns as
#' `x_train`.
#'
#' @param train_df Training data frame.
#' @param new_df   New (test/validation) data frame.
#' @param feature_cols Character vector of feature column names.
#'
#' @return A list with `x_train` and `x_new`, both numeric matrices.
#' @export
remake_X <- function(train_df, new_df, feature_cols) {
  tr <- train_df[, feature_cols, drop = FALSE]
  nw <- new_df[,  feature_cols, drop = FALSE]
  
  fcols <- names(tr)[vapply(tr, is.factor, logical(1L))]
  for (cn in fcols)
    nw[[cn]] <- factor(nw[[cn]], levels = levels(tr[[cn]]))
  
  tt    <- stats::terms(~ ., data = tr)
  mf_tr <- stats::model.frame(tt, data = tr, na.action = stats::na.pass)
  X_tr  <- stats::model.matrix(tt, data = mf_tr)[, -1L, drop = FALSE]
  
  tt    <- stats::terms(~ ., data = nw)
  mf_nw <- stats::model.frame(tt, data = nw, na.action = stats::na.pass)
  X_nw  <- stats::model.matrix(tt, data = mf_nw)[, -1L, drop = FALSE]
  
  colnames(X_tr) <- make.names(colnames(X_tr), unique = TRUE)
  colnames(X_nw) <- make.names(colnames(X_nw), unique = TRUE)
  
  miss <- setdiff(colnames(X_tr), colnames(X_nw))
  if (length(miss))
    X_nw <- cbind(X_nw,
                  matrix(0, nrow(X_nw), length(miss),
                         dimnames = list(NULL, miss)))
  X_nw <- X_nw[, colnames(X_tr), drop = FALSE]
  
  keep_var <- apply(X_tr, 2L, function(z) stats::sd(z) > 0)
  if (any(!keep_var)) {
    X_tr <- X_tr[, keep_var, drop = FALSE]
    X_nw <- X_nw[, keep_var, drop = FALSE]
  }
  
  if (ncol(X_tr) > 0L) {
    q   <- qr(X_tr)
    piv <- q$pivot[seq_len(q$rank)]
    X_tr <- X_tr[, piv, drop = FALSE]
    X_nw <- X_nw[, colnames(X_tr), drop = FALSE]
  }
  
  storage.mode(X_tr) <- "double"
  storage.mode(X_nw) <- "double"
  
  list(x_train = X_tr, x_new = X_nw)
}


#' Build the output path for storing CIF predictions
#'
#' @param out_dir Root output directory.
#' @param model  Model key string.
#' Stratified k-fold CV split object
#'
#' Creates an `rsample` v-fold cross-validation object stratified by the
#' event column so that each fold contains a representative mix of event
#' types.
#'
#' @param df A data frame containing an `event` column.
#' @param time  (Unused; kept for signature compatibility.)
#' @param event (Unused; `df$event` is used directly.)
#' @param folds Number of folds (default 5).
#'
#' @return An `rsample::vfold_cv` object.
#' @export
stratification <- function(df, time, event, folds = 5) {
  df$event <- as.factor(df$event)
  rsample::vfold_cv(df, v = folds, strata = "event")
}


#' @title Build storage path for CIF outputs
#' @description Creates the directory `out_dir/<model>/outer_<fold>/` and
#'   returns the path.
#'
#' @param out_dir  Root output directory.
#' @param model    Model key string.
#' @param fold     Outer fold index.
#'
#' @return The path (created on disk if needed).
#' @export
build_store_paths_r <- function(out_dir, model, fold) {
  p <- file.path(out_dir, model, sprintf("outer_%d", fold))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  p
}


#' Save CIF predictions to a Parquet file
#'
#' Flattens the 3-D CIF array into a long data frame with columns
#' `row_id`, `cause`, `time`, and `cif`, then writes it as a Parquet file.
#'
#' @param out_dir Directory in which `cif.parquet` is written.
#' @param cif     3-D array with dimensions `[n, K, Tm]`.
#' @param times   Numeric vector of evaluation times (length `Tm`).
#' @param row_ids Character vector of subject identifiers (length `n`).
#' @param causes  Integer vector of cause codes (length `K`).
#'
#' @return Invisibly `NULL`.
#' @export
save_cif_r <- function(out_dir, cif, times, row_ids, causes) {
  n  <- dim(cif)[1L]; K <- dim(cif)[2L]; Tm <- dim(cif)[3L]
  df <- data.frame(
    row_id = rep(row_ids, each = K * Tm),
    cause  = rep(rep(causes, each = Tm), times = n),
    time   = rep(times, times = n * K),
    cif    = as.vector(cif)
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  arrow::write_parquet(df, file.path(out_dir, "cif.parquet"))
  invisible(NULL)
}


#' Extract penalised Fine-Gray regression coefficients
#'
#' Convenience helper to pull the named coefficient vector from a fitted
#' FGRP model object (as returned by the registered `"FGRP"` model).
#'
#' @param fit_obj A `cr_model_fgrp` object returned by the FGRP `fit()`
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


#' Simulate a competing risks dataset
#'
#' Generates synthetic data following a DeepHit-like data-generating
#' process with two competing causes and configurable covariate blocks.
#'
#' @param n         Number of subjects (default 30 000).
#' @param p_block   Number of covariates per linear predictor block
#'   (default 4).
#' @param gammas    Shared coefficient magnitude for all blocks (default 10).
#' @param cens_frac Fraction of subjects to be censored (default 0.5).
#' @param cap       Optional symmetric cap applied to the linear predictor to
#'   avoid numerical overflow.
#' @param seed      Random seed (default 123).
#'

#' Trapezoidal numerical integration
#'
#' @param x Numeric vector of x-values (must be sorted).
#' @param y Numeric vector of y-values (same length as `x`).
#'
#' @return A single numeric value.
#' @export
trapezoidal.integration <- function(x, y) {
  if (length(x) != length(y))
    stop("`x` and `y` must have the same length.")
  sum(diff(x) * (utils::head(y, -1L) + utils::tail(y, -1L)) / 2)
}


#' Compute the restricted mean lifetime (RMLT) from a CIF array
#'
#' Integrates each subject's CIF up to `tau` using trapezoidal integration.
#'
#' @param out   3-D CIF array with dimensions `[n, K, Tm]`.
#' @param times Numeric vector of evaluation times (length `Tm`).
#' @param tau   Upper integration limit (default `max(times)`).
#'
#' @return An `[n x K]` numeric matrix of per-subject RMLTs.
#' @export
compute_rmlt_from_cif <- function(out, times, tau = max(times)) {
  d     <- dim(out)
  n     <- d[1L]; K <- d[2L]
  idx_t <- which(times <= tau)
  
  if (length(idx_t) < 2L) {
    res <- matrix(0, nrow = n, ncol = K)
    colnames(res) <- paste0("cause_", seq_len(K))
    return(res)
  }
  
  times_use <- times[idx_t]
  res <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    cif_mat    <- matrix(out[, k, idx_t, drop = FALSE], nrow = n)
    res[, k]   <- apply(cif_mat, 1L,
                        trapezoidal.integration, x = times_use)
  }
  colnames(res) <- paste0("cause_", seq_len(K))
  res
}


#' Build the imputation cache for a benchmark directory
#'
#' Pre-computes and saves imputed train/test Parquet files for all outer
#' and inner folds.  The resulting directory can be used as a drop-in
#' replacement for the original benchmark directory.
#'
#' @param out_dir        Path to the benchmark output directory.
#' @param impute_method  One of `"mice"` or `"median_mode"`.
#' @param seed           Integer random seed.
#' @param m              Number of imputed datasets (MICE only).
#' @param maxit          MICE iterations.
#' @param drop_cols_after_impute Optional columns to drop after imputation.
#' @param overwrite      If `TRUE`, overwrite existing Parquet files.
#' @param verbose        Print progress messages.
#'
#' @return Invisibly, the path to the cache directory.
#' @export
build_imputation_cache <- function(out_dir                = "../BenchResults",
                                   impute_method          = c("mice", "median_mode"),
                                   seed                   = 123L,
                                   m                      = 1L,
                                   maxit                  = 20L,
                                   drop_cols_after_impute = NULL,
                                   overwrite              = FALSE,
                                   verbose                = TRUE) {
  impute_method <- match.arg(impute_method)
  
  manifest_path <- file.path(out_dir, "manifest.json")
  schema_path   <- file.path(out_dir, "data_types.json")
  times_path    <- file.path(out_dir, "times.json")
  
  man    <- jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
  schema <- jsonlite::fromJSON(schema_path,   simplifyVector = TRUE)
  
  outer_folds  <- as.integer(man$outer_folds)
  feature_cols <- setdiff(man$features, "row_id")
  
  cache_dir <- file.path(
    dirname(out_dir),
    paste0(basename(out_dir), "_imputed_", impute_method)
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  file.copy(manifest_path, file.path(cache_dir, "manifest.json"),
            overwrite = TRUE)
  if (file.exists(schema_path))
    file.copy(schema_path, file.path(cache_dir, "data_types.json"),
              overwrite = TRUE)
  if (file.exists(times_path))
    file.copy(times_path, file.path(cache_dir, "times.json"),
              overwrite = TRUE)
  
  add_row_id <- function(df)
    data.frame(row_id = rownames(df), df, check.names = FALSE)
  
  write_if_needed <- function(df, path) {
    if (!overwrite && file.exists(path)) return(invisible(FALSE))
    arrow::write_parquet(df, path, compression = "zstd")
    invisible(TRUE)
  }
  
  write_imp_pair <- function(imp, out_dir_split, test_prefix = "test") {
    dir.create(out_dir_split, recursive = TRUE, showWarnings = FALSE)
    for (k in seq_len(imp$m)) {
      tr_k <- add_row_id(imp$imputations[[k]]$train)
      te_k <- add_row_id(imp$imputations[[k]]$test)
      write_if_needed(tr_k,
                      file.path(out_dir_split,
                                sprintf("train_imp_%d.parquet", k)))
      write_if_needed(te_k,
                      file.path(out_dir_split,
                                sprintf("%s_imp_%d.parquet", test_prefix, k)))
    }
  }
  
  for (v in seq_len(outer_folds)) {
    fold_dir <- file.path(out_dir, sprintf("outer_%d", v))
    train <- apply_data_types(
      as.data.frame(arrow::read_parquet(file.path(fold_dir, "train.parquet"))),
      schema = schema
    )
    test <- apply_data_types(
      as.data.frame(arrow::read_parquet(file.path(fold_dir, "test.parquet"))),
      schema = schema
    )
    rownames(train) <- as.character(train$row_id); train$row_id <- NULL
    rownames(test)  <- as.character(test$row_id);  test$row_id  <- NULL
    
    out_fold_cache <- file.path(cache_dir, sprintf("outer_%d", v))
    if (verbose) message("[impute->parquet] ", impute_method, " outer ", v)
    imp_outer <- impute_split(
      train = train, test = test, feature_cols = feature_cols,
      drop_cols_after_impute = drop_cols_after_impute,
      method = impute_method, m = m, maxit = maxit, seed = seed + 1000L + v
    )
    write_imp_pair(imp_outer, out_fold_cache, test_prefix = "test")
    
    inner_dir  <- file.path(fold_dir, "inner")
    if (!dir.exists(inner_dir)) next
    
    inner_dirs <- list.dirs(inner_dir, recursive = FALSE, full.names = TRUE)
    inner_dirs <- inner_dirs[grepl("fold_[0-9]+$", basename(inner_dirs))]
    
    for (j in seq_along(inner_dirs)) {
      idx_path <- file.path(inner_dirs[j], "indices.json")
      idx      <- jsonlite::fromJSON(idx_path, simplifyVector = TRUE)
      tr_ids   <- as.character(idx$train_row_id)
      va_ids   <- as.character(idx$val_row_id)
      
      train_in <- train[tr_ids, , drop = FALSE]
      val_in   <- train[va_ids, , drop = FALSE]
      
      in_cache_dir <- file.path(out_fold_cache, "inner",
                                basename(inner_dirs[j]))
      dir.create(in_cache_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(idx_path, file.path(in_cache_dir, "indices.json"),
                overwrite = TRUE)
      
      if (verbose)
        message("[impute->parquet] ", impute_method,
                " outer ", v, " inner ", j)
      imp_inner <- impute_split(
        train = train_in, test = val_in, feature_cols = feature_cols,
        drop_cols_after_impute = drop_cols_after_impute,
        method = impute_method, m = m, maxit = maxit,
        seed = seed + 2000L + v * 100L + j
      )
      write_imp_pair(imp_inner, in_cache_dir, test_prefix = "val")
    }
  }
  
  invisible(cache_dir)
}


#' Summarise nested cross-validation results across outer folds
#'
#' @param results A list of per-fold result lists as returned by
#'   [nested_cv_from_bench()].
#' @param num_causes Integer; number of competing causes.
#' @param times Numeric vector of evaluation time points.
#'
#' @return A data frame with one row per model × fold × cause × time,
#'   containing Brier score, IBS, AUC, C-index (pec), C-index (survMetrics),
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
        bs         <- get_or_na(res_fold$cause_bs[[cause]],          length(times))
        ibs        <- get_or_na(res_fold$cause_ibs[[cause]],         length(times))
        auc        <- get_or_na(res_fold$cause_auc[[cause]],         length(times))
        cidx_pec   <- get_or_na(res_fold$cause_cindex_pec[[cause]],  length(times))
        cidx_survM <- if (!is.null(res_fold$cause_cindex_survM[[cause]]))
          res_fold$cause_cindex_survM[[cause]]
        else
          list(rep(NA_real_, length(times)), rep(NA_real_, length(times)))
        calib_meas <- if (!is.null(res_fold$cause_calib_measures[[cause]]))
          res_fold$cause_calib_measures[[cause]]
        else
          list(ICI = NA_real_, E50 = NA_real_, E90 = NA_real_,
               Emax = NA_real_, RSB = NA_real_)
        
        tab[[length(tab) + 1L]] <- data.frame(
          model      = model,
          fold       = v,
          cause      = cause,
          bs         = bs,
          ibs        = ibs,
          auc        = auc,
          cidx_pec   = cidx_pec,
          cidx_survM = cidx_survM[[1L]],
          survM_eval = cidx_survM[[2L]],
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
  if ("cidx_survM" %in% metrics)
    out[["survM_eval"]] <- agg_one("survM_eval")
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
                                            "cidx_pec", "cidx_survM",
                                            "ICI", "E50", "E90",
                                            "Emax", "RSB"),
                               cause    = NULL,
                               plot_tau = TRUE) {
  for (metric in metrics) {
    ylimit <- c(0, 1)
    df     <- agg[[metric]]
    if (!is.null(cause)) df <- df[df$cause == cause, , drop = FALSE]
    
    tau <- tau_min <- tau_max <- tau_med <- tau_mu <- NULL
    
    if (metric == "cidx_survM" && !is.null(agg$survM_eval)) {
      tau_df  <- agg$survM_eval
      tau_min <- min(tau_df$x[, "lwr"],    na.rm = TRUE)
      tau_max <- max(tau_df$x[, "upr"],    na.rm = TRUE)
      tau_med <- unique(tau_df$x[, "median"])
      tau_mu  <- unique(tau_df$x[, "mean"])
    } else if (metric %in% c("bs", "auc", "cidx_pec",
                             "ICI", "E50", "E90", "Emax", "RSB")) {
      tau <- df$time
      if (length(unique(tau)) == 1L && metric == "cidx_pec") {
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
    
    if (plot_tau && metric == "cidx_survM" &&
        !is.null(tau_min) && !is.null(tau_max)) {
      p <- p + ggplot2::annotate("rect",
                                 xmin = tau_min, xmax = tau_max,
                                 ymin = -Inf, ymax = Inf,
                                 fill = "grey90", alpha = 0.05)
    }
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
  cfg_rows       <- list()
  ibs_rows       <- list()
  ibs_time_rows  <- list()
  
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
      cfg_rows[[length(cfg_rows) + 1L]] <- cfg_df
      
      ibs_vec <- as.numeric(r$metrics$ibs_scores[[ci]])
      ibs_rows[[length(ibs_rows) + 1L]] <- data.frame(
        fold = fold, model = model, cause = k,
        ibs  = utils::tail(ibs_vec, 1L),
        stringsAsFactors = FALSE
      )
      
      Tm     <- length(ibs_vec)
      t_used <- if (length(times) == Tm) times else seq_len(Tm)
      ibs_time_rows[[length(ibs_time_rows) + 1L]] <- data.frame(
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
  
  ibs_summary <- if (nrow(ibs_per_fold_cause) > 0L) {
    key <- interaction(ibs_per_fold_cause$model,
                       ibs_per_fold_cause$cause, drop = TRUE)
    spl <- split(ibs_per_fold_cause$ibs, key)
    stats <- lapply(spl, function(x)
      c(mean = mean(x), sd = stats::sd(x), median = stats::median(x)))
    out <- do.call(rbind, stats)
    keys   <- names(spl)
    parts  <- read.table(text = keys, sep = ".",
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
        out <- aperm(array(df_cif$cif, dim = c(K, Tm, n)), c(3L, 1L, 2L))
      } else {
        out <- array(df_cif$cif, dim = c(n, K, Tm))
      }
      
      if (sample_curves) {
        set.seed(seed)
        ids <- sample.int(n, min(20L, n))
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
      
      res[[model]][[v]] <- score_from_cifs(
        out        = out,
        test       = df_test,
        times      = times,
        causes     = causes,
        time_col   = "time",
        status_col = "event",
        metrics    = c("brier", "auc", "calib_measures"),
        summary    = "ibs",
        cens.method = "ipcw",
        cens.model  = "km"
      )
      
      all_cifs[[model]][[v]] <- out
      all_test[[v]]          <- df_test
      
      all_rmlt[[model]][[v]] <- compute_rmlt_from_cif(out, times,
                                                      tau = max(times))
      res_rmlt[[model]][[v]] <- score_from_rmlt(
        rmlt       = all_rmlt[[model]][[v]],
        test       = df_test,
        times      = max(times),
        causes     = causes,
        time_col   = "time",
        status_col = "event",
        metrics    = "cidx_pec"
      )
    }
    
    all_cifs[[model]] <- abind::abind(all_cifs[[model]], along = 1L)
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
      cens.code        = 0L,
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