#' @import S7
#' @importFrom Rcpp sourceCpp
#' @useDynLib gPIPE, .registration = TRUE
NULL

.onLoad <- function(...) {
  S7::S7_on_load()
}

.onUnload <- function(...) {
  S7::S7_on_unload()
}

S7::S7_on_build()
