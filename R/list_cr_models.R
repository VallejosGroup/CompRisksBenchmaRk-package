#' List all registered model keys
#'
#' @return A character vector of registered model keys (may be empty).
#' @export
#'
#' @examples
#' list_cr_models()
list_cr_models <- function() ls(envir = .cr_models)
