#' @title CIF Array and Design Matrix Utilities
#' @description Functions for pooling, saving, and storing CIF prediction
#'   arrays, and for building aligned numeric design matrices.
#' @name cif_io
NULL


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
  if (m > 1)
    for (i in 2:m) pooled <- pooled + cif_list[[i]]
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

  fcols <- names(tr)[vapply(tr, is.factor, logical(1))]
  for (cn in fcols)
    nw[[cn]] <- factor(nw[[cn]], levels = levels(tr[[cn]]))

  tt    <- stats::terms(~ ., data = tr)
  mf_tr <- stats::model.frame(tt, data = tr, na.action = stats::na.pass)
  X_tr  <- stats::model.matrix(tt, data = mf_tr)[, -1, drop = FALSE]

  tt    <- stats::terms(~ ., data = nw)
  mf_nw <- stats::model.frame(tt, data = nw, na.action = stats::na.pass)
  X_nw  <- stats::model.matrix(tt, data = mf_nw)[, -1, drop = FALSE]

  colnames(X_tr) <- make.names(colnames(X_tr), unique = TRUE)
  colnames(X_nw) <- make.names(colnames(X_nw), unique = TRUE)

  miss <- setdiff(colnames(X_tr), colnames(X_nw))
  if (length(miss))
    X_nw <- cbind(X_nw,
                  matrix(0, nrow(X_nw), length(miss),
                         dimnames = list(NULL, miss)))
  X_nw <- X_nw[, colnames(X_tr), drop = FALSE]

  keep_var <- apply(X_tr, 2, function(z) stats::sd(z) > 0)
  if (any(!keep_var)) {
    X_tr <- X_tr[, keep_var, drop = FALSE]
    X_nw <- X_nw[, keep_var, drop = FALSE]
  }

  if (ncol(X_tr) > 0) {
    q   <- qr(X_tr)
    piv <- q$pivot[seq_len(q$rank)]
    X_tr <- X_tr[, piv, drop = FALSE]
    X_nw <- X_nw[, colnames(X_tr), drop = FALSE]
  }

  storage.mode(X_tr) <- "double"
  storage.mode(X_nw) <- "double"

  list(x_train = X_tr, x_new = X_nw)
}


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
#' @param out_dir Directory in which `cif.parquet` is written.
#' @param cif     3-D array with dimensions `[n, K, Tm]`.
#' @param times   Numeric vector of evaluation times (length `Tm`).
#' @param row_ids Character vector of subject identifiers (length `n`).
#' @param causes  Integer vector of cause codes (length `K`).
#'
#' @return Invisibly `NULL`.
#' @export
save_cif_r <- function(out_dir, cif, times, row_ids, causes) {
  n  <- dim(cif)[1]; K <- dim(cif)[2]; Tm <- dim(cif)[3]
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
