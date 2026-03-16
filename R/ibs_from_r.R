#' @title Performance Metrics for Competing Risks Models
#' @description Functions for computing performance metrics for competing risks
#'   survival models.
#' @name ibs_from_r
NULL


#' IBS wrapper for use from Python or external processes
#'
#' Computes the Integrated Brier Score for each cause in a 3-D CIF array
#' using `riskRegression::Score()`.
#'
#' @param out         3-D numeric array `[n, K, Tm]`.
#' @param e_va        Integer vector of observed event codes (length `n`).
#' @param t_va        Numeric vector of observed times (length `n`).
#' @param times       Numeric vector of evaluation times (length `Tm`).
#' @param causes      Integer vector of cause codes (length `K`).
#' @param time_col    Column name for time (default `"time"`).
#' @param status_col  Column name for status (default `"event"`).
#' @param cens.method Passed to `Score()` (default `"ipcw"`).
#' @param cens.model  Passed to `Score()` (default `"km"`).
#' @param cens.code   Censoring code (default `0`).
#' @param force_sequential If `TRUE`, force sequential execution and
#'   single-threaded BLAS (default `TRUE`).
#' @param quiet       Suppress output (default `TRUE`).
#' @param log_file    Optional path to a log file.
#' @param debug       Return debug information? (default `FALSE`).
#' @param debug_n     Number of rows to include in debug output.
#'
#' @return A numeric matrix `[K, Tm]` of IBS values, or a debug list if
#'   `debug = TRUE`.
#' @export
ibs_from_r <- function(out, e_va, t_va, times, causes,
                        time_col          = "time",
                        status_col        = "event",
                        cens.method       = "ipcw",
                        cens.model        = "km",
                        cens.code         = 0L,
                        force_sequential  = TRUE,
                        quiet             = TRUE,
                        log_file          = NULL,
                        debug             = FALSE,
                        debug_n           = 5L) {
  if (isTRUE(force_sequential)) {
    foreach::registerDoSEQ()
    data.table::setDTthreads(1L)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
  }

  log <- function(...) {
    msg <- paste0(paste(..., collapse = ""), "\n")
    if (!is.null(log_file)) cat(msg, file = log_file, append = TRUE)
    else if (!isTRUE(quiet)) cat(msg)
  }

  d  <- dim(out)
  if (length(d) != 3L) stop("`out` must be a 3D array (n, K, T).")
  n  <- d[1L]; K <- d[2L]; Tm <- d[3L]
  log("dim(out) = ", paste(d, collapse = " x "))
  if (length(times)  != Tm) stop("length(times) must equal dim(out)[3].")
  if (length(causes) != K)  stop("length(causes) must equal dim(out)[2].")

  times <- as.numeric(times)
  ord   <- order(times)
  if (any(ord != seq_along(times))) {
    out   <- out[, , ord, drop = FALSE]
    times <- times[ord]
  }
  if (any(duplicated(times))) stop("`times` must be unique.")

  test        <- data.frame(as.numeric(t_va), as.integer(e_va))
  names(test) <- c(time_col, status_col)
  f           <- stats::as.formula(
    sprintf("prodlim::Hist(%s, %s) ~ 1", time_col, status_col)
  )

  run_quiet <- function(expr) {
    if (!isTRUE(quiet)) return(eval.parent(substitute(expr)))
    res <- NULL
    invisible(utils::capture.output({
      res <- withCallingHandlers(
        eval.parent(substitute(expr)),
        warning = function(w) invokeRestart("muffleWarning")
      )
    }))
    res
  }

  ibs_mat <- matrix(NA_real_, nrow = length(causes), ncol = length(times),
                    dimnames = list(as.character(causes),
                                   as.character(times)))

  for (i in seq_len(K)) {
    k      <- as.integer(causes[i])
    M      <- out[, i, , drop = FALSE]
    dim(M) <- c(n, Tm)
    M      <- as.matrix(M)
    storage.mode(M) <- "double"
    colnames(M) <- format(times, scientific = FALSE, trim = TRUE)

    k_name     <- paste0("cif_cause_", k)
    preds_list <- stats::setNames(list(M), k_name)

    sc <- run_quiet(
      riskRegression::Score(
        object       = preds_list,
        formula      = f,
        data         = test,
        metrics      = "brier",
        summary      = "ibs",
        times        = times,
        cause        = k,
        cens.method  = cens.method,
        cens.model   = cens.model,
        cens.code    = cens.code,
        split.method = "none",
        se.fit       = FALSE,
        conf.int     = FALSE,
        progress.bar = FALSE,
        verbose      = 0,
        censoring.save.memory = TRUE,
        keep         = NULL,
        plots        = NULL
      )
    )

    ibs_val        <- sc$Brier$score[sc$Brier$score$model == k_name, "IBS"]
    ibs_mat[i, ]   <- as.numeric(ibs_val)
    rm(sc, M, preds_list); gc(FALSE)
  }

  if (isTRUE(debug)) {
    out_slice <- out[seq_len(min(n, debug_n)),
                     seq_len(min(K, debug_n)),
                     seq_len(min(Tm, debug_n)), drop = FALSE]
    return(list(
      ibs   = ibs_mat,
      debug = list(dim_out   = d, dim_slice = dim(out_slice),
                   times_head = utils::head(times, debug_n),
                   causes = causes, test_head = utils::head(test, debug_n),
                   event_table = table(test[[status_col]], useNA = "ifany"))
    ))
  }
  ibs_mat
}

