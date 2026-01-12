#' @keywords internal
"_PACKAGE"

# Quiet R CMD check notes about global variables
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("w", "ipw_weights", "weights"))
}

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib ordinal2phase, .registration = TRUE
## usethis namespace: end
NULL
