#' @title CIF Array and Storage Utilities
#' @description Functions for saving and storing CIF prediction arrays.
#' @name cif_io
NULL



#' Build storage path for CIF outputs
#'
#' Creates the directory `out_dir/<model>/outer_<fold>/` and returns the path.
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
#' @param out_dir       Directory in which `cif.parquet` is written.
#' @param cif           3-D array with dimensions `[n, K, Tm]`.
#' @param cif_time_grid Numeric vector of evaluation times (length `Tm`).
#' @param row_ids       Character vector of subject identifiers (length `n`).
#' @param causes        Integer vector of cause codes (length `K`).
#'
#' @return Invisibly `NULL`.
#' @export
save_cif_r <- function(out_dir, cif, cif_time_grid, row_ids, causes) {
  n  <- dim(cif)[1]; K <- dim(cif)[2]; Tm <- dim(cif)[3]
  df <- data.frame(
    row_id = rep(row_ids, each = K * Tm),
    cause  = rep(rep(causes, each = Tm), times = n),
    time   = rep(cif_time_grid, times = n * K),
    cif    = as.vector(cif)
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  arrow::write_parquet(df, file.path(out_dir, "cif.parquet"))
  invisible(NULL)
}

#' Evaluate metrics from stored nested CV predictions
#'
#' Reads the CIF predictions written by [nested_cv_from_bench()] and the
#' original test data from the benchmark directory, then calls
#' [compute_metrics()] for each outer fold.
#'
#' @param out_dir       Path to the benchmark directory used in
#'   [nested_cv_from_bench()].
#' @param model_key     Character key of the model whose predictions to
#'   evaluate.
#' @param pred_horizons Numeric vector of prediction horizons passed to
#'   [compute_metrics()].
#' @param metrics       Character vector of metrics to compute; passed to
#'   [compute_metrics()].
#' @param collapse_as_df Logical; passed to [compute_metrics()].
#'
#' @return A list of length `outer_folds`, each element being the output of
#'   [compute_metrics()] for that fold.
#' @export
eval_nested_cv <- function(out_dir,
                           model_key,
                           pred_horizons,
                           metrics        = c("Brier", "IBS"),
                           collapse_as_df = TRUE) {

  # Load CV metadata (includes time_var, event_var, id_var, outer_folds, etc.)
  cv_meta <- jsonlite::fromJSON(
    file.path(out_dir, "cv_splits_metadata.json"),
    simplifyVector = TRUE
  )

  outer_folds <- as.integer(cv_meta$outer_folds)
  time_var    <- cv_meta$time_var
  event_var   <- cv_meta$event_var
  id_var      <- cv_meta$id_var
  causes      <- as.integer(cv_meta$causes)

  results <- vector("list", outer_folds)

  for (v in seq_len(outer_folds)) {
    fold_dir <- file.path(out_dir, sprintf("outer_%d", v))
    cif_path <- file.path(out_dir, model_key, sprintf("outer_%d", v), "cif.parquet")

    if (!file.exists(cif_path))
      stop(sprintf("CIF file not found for fold %d: '%s'", v, cif_path), call. = FALSE)

    # Load test data and reconstruct cr_data object
    test_df <- .apply_var_types(
      as.data.frame(arrow::read_parquet(file.path(fold_dir, "test.parquet"))),
      cv_meta$var_types,
      cv_meta$var_names
    )
    cr_test <- cr_data(test_df, time_var = time_var, event_var = event_var,
                       id_var = id_var)

    # Read stored CIF predictions and reconstruct 3-D array [n, K, Tm]
    cif_df    <- as.data.frame(arrow::read_parquet(cif_path))
    row_ids   <- unique(cif_df$row_id)
    cif_times <- sort(unique(cif_df$time))
    n  <- length(row_ids)
    K  <- length(causes)
    Tm <- length(cif_times)

    cif_arr <- array(NA_real_, dim = c(n, K, Tm))
    for (ci in seq_along(causes)) {
      sub <- cif_df[cif_df$cause == causes[ci], ]
      sub <- sub[order(match(sub$row_id, row_ids), match(sub$time, cif_times)), ]
      cif_arr[, ci, ] <- matrix(sub$cif, nrow = n, ncol = Tm, byrow = FALSE)
    }

    results[[v]] <- compute_metrics(
      cr_test,
      cif           = list(cif = cif_arr, time_grid = cif_times),
      pred_horizons = pred_horizons,
      metrics       = metrics,
      collapse_as_df = collapse_as_df
    )
  }

  results
}
