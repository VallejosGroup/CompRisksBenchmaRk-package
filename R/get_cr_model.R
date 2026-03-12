#' Retrieve a registered competing risks model
#'
#' @param key Character string matching a previously registered model key.
#'
#' @return A named list with elements `fit`, `predict_cif`, and `info`.
#' @export
#'
#' @examples
#' \dontrun{
#' mdl <- get_cr_model("FGR")
#' mdl$info()
#' }
get_cr_model <- function(key) {
  m <- get0(key, envir = .cr_models, inherits = FALSE)
  if (is.null(m)) stop("Model not registered: ", key, call. = FALSE)
  m
}
