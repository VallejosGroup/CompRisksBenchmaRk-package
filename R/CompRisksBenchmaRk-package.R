#' CompRisksBenchmaRk: A Benchmarking Framework for Competing Risks Models
#'
#' @description
#' \pkg{CompRisksBenchmaRk} provides a unified workflow for fitting,
#' predicting, and benchmarking competing risks survival models.  It
#' includes:
#'
#' \describe{
#'   \item{Model registry}{Register any custom competing risks model with
#'     [register_cr_model()], retrieve it with [get_cr_model()], and list
#'     all registered models with [list_cr_models()].}
#'   \item{Built-in models}{Fine-Gray (`"FGR"`), penalised Fine-Gray
#'     (`"FGRP"`), cause-specific Cox PH (`"csCPH"`), and random survival
#'     forest (`"RSF"`) are registered automatically on package load.}
#'   \item{Nested cross-validation}{[nested_cv_from_bench()] runs a full
#'     nested CV pipeline from pre-built splits created by
#'     [create_nested_splits()].}
#'   \item{Performance metrics}{Time-varying and integrated Brier scores,
#'     AUC, C-index (pec and SurvMetrics), calibration curves, and
#'     observed/expected ratios.}
#'   \item{Clinical utility}{Decision curve analysis via
#'     [ClinicalUtility()].}
#'   \item{Imputation}{Median/mode and MICE-based missing-data handling
#'     fully integrated into the CV pipeline.}
#' }
#'
#' @keywords internal
"_PACKAGE"
