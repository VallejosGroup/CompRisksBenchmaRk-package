#' @title Internal data loading utilities for nested CV
#' @description Internal helpers for reading Parquet fold data and inner split
#'   indices from the benchmark directory layout produced by
#'   [create_nested_splits()].
#' @name utils_load
NULL


#' @noRd
.apply_var_types <- function(df, var_types, var_names) {
  for (i in seq_along(var_names)) {
    cn <- var_names[[i]]
    if (!cn %in% names(df)) next
    df[[cn]] <- switch(var_types[[i]],
                       numeric   = as.numeric(df[[cn]]),
                       integer   = as.integer(df[[cn]]),
                       logical   = as.logical(df[[cn]]),
                       factor    = as.factor(df[[cn]]),
                       as.character(df[[cn]])
    )
  }
  df
}


#' @noRd
.detect_m_from_dir <- function(dir) {
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
.load_outer_data <- function(fold_dir, metadata) {
  if (file.exists(file.path(fold_dir, "train_imp_1.parquet")) &&
      file.exists(file.path(fold_dir, "test_imp_1.parquet"))) {
    mm  <- .detect_m_from_dir(fold_dir)
    out <- vector("list", mm)
    for (k in seq_len(mm)) {
      trp <- file.path(fold_dir, sprintf("train_imp_%d.parquet", k))
      tep <- file.path(fold_dir, sprintf("test_imp_%d.parquet",  k))
      out[[k]] <- list(
        train = .apply_var_types(as.data.frame(arrow::read_parquet(trp)), metadata$var_types, metadata$var_names),
        test  = .apply_var_types(as.data.frame(arrow::read_parquet(tep)), metadata$var_types, metadata$var_names))
    }
    return(out)
  }
  tr_path <- file.path(fold_dir, "train.parquet")
  te_path <- file.path(fold_dir, "test.parquet")
  list(list(
    train = .apply_var_types(as.data.frame(arrow::read_parquet(tr_path)), metadata$var_types, metadata$var_names),
    test  = .apply_var_types(as.data.frame(arrow::read_parquet(te_path)), metadata$var_types, metadata$var_names)))
}


#' @noRd
.load_inner_data <- function(fold_dir, metadata, outer_train_df) {
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
      mm <- .detect_m_from_dir(idir)
      jj <- vector("list", mm)
      for (k in seq_len(mm)) {
        trp <- file.path(idir, sprintf("train_imp_%d.parquet", k))
        vap <- file.path(idir, sprintf("val_imp_%d.parquet",   k))
        jj[[k]] <- list(
          train = .apply_var_types(as.data.frame(arrow::read_parquet(trp)), metadata$var_types, metadata$var_names),
          val   = .apply_var_types(as.data.frame(arrow::read_parquet(vap)), metadata$var_types, metadata$var_names))
      }
      inner_data[[j]] <- jj
      next
    }
    idx_path <- file.path(idir, "indices.json")
    idx      <- jsonlite::fromJSON(idx_path, simplifyVector = TRUE)
    tr_ids   <- intersect(as.character(idx$train_row_id), rownames(outer_train_df))
    va_ids   <- intersect(as.character(idx$val_row_id),   rownames(outer_train_df))
    inner_data[[j]] <- list(list(
      train = outer_train_df[tr_ids, , drop = FALSE],
      val   = outer_train_df[va_ids, , drop = FALSE]
    ))
  }
  
  list(inner_dirs = inner_dirs, inner_data = inner_data)
}