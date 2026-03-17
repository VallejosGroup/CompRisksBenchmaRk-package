#' @title Internal data loading utilities for nested CV
#' @description Internal helpers for reading Parquet fold data and inner split
#'   indices from the benchmark directory layout produced by
#'   [create_nested_splits()].
#' @name utils_load
NULL


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
