#' @title Imputation Utilities
#' @description Functions for detecting and filling missing values in
#'   train/test splits, supporting median/mode and MICE strategies.
#' @name imputation
NULL


#' @noRd
split_detect_NA <- function(train, test, cols) {
  anyNA(train[, cols, drop = FALSE]) ||
    anyNA(test[, cols, drop = FALSE])
}


#' @noRd
mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA)
  tab <- table(x)
  names(tab)[which.max(tab)]
}


#' @noRd
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


#' @noRd
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


#' @noRd
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

  keep_cols   <- setdiff(names(train), feature_cols)
  imputations <- vector("list", m)

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
                          m      = 1,
                          maxit  = 100,
                          seed   = 1,
                          drop_cols_after_impute = NULL) {
  method <- match.arg(method)

  if (method == "mice") {
    cat("Default number of datasets m: ", m, ", and iterations: ", maxit,
        " seed: ", seed, "\n")
  }

  if (!split_detect_NA(train, test, feature_cols) || method == "none") {
    return(list(m = 1, imputations = list(list(train = train, test = test))))
  }

  if (method == "median_mode") {
    imp <- impute_median_mode_single(
      train = train, test = test,
      feature_cols = feature_cols,
      drop_cols_after_impute = drop_cols_after_impute
    )
    return(list(m = 1,
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
                                    seed                   = 123,
                                    m                      = 1,
                                    maxit                  = 20,
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
      method = impute_method, m = m, maxit = maxit, seed = seed + 1000 + v
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
        seed = seed + 2000 + v * 100 + j
      )
      write_imp_pair(imp_inner, in_cache_dir, test_prefix = "val")
    }
  }

  invisible(cache_dir)
}
