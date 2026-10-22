#' Generate All Monotonic Binary Matrices
#'
#' This function generates all possible binary (0/1) matrices of a given dimension
#' that satisfy monotonicity constraints in both row and column directions.
#' These matrices can be used to represent valid dose configurations in dose-finding studies
#' where toxicity is assumed to be monotonic with respect to dose levels.
#'
#' @param I Integer specifying the number of rows
#' @param J Integer specifying the number of columns
#' @param direction Character string specifying the direction of monotonicity:
#'    \itemize{
#'      \item "increasing": Values increase or stay the same as row/column indices increase
#'      \item "decreasing": Values decrease or stay the same as row/column indices increase
#'    }
#'
#' @return A list of all possible monotonic binary matrices of dimension I×J
#'
#' @details
#' For dose-finding studies, these matrices can represent valid dose-toxicity relationships:
#'   \itemize{
#'     \item 1 represents acceptable toxicity
#'     \item 0 represents unacceptable toxicity
#'   }
#' 
#' In the "increasing" direction, a dose is acceptable only if all lower doses are also acceptable.
#' In the "decreasing" direction, a dose is unacceptable only if all higher doses are also unacceptable.
#'
#' The matrices come from the up-set enumeration of a [DoseGrid], in the
#' order of the original recursive version (cells read row by row, 0 before 1).
#'
#' @examples
#' # Generate all 2x2 matrices with increasing monotonicity
#' matrices_inc <- monotonic_matrices(2, 2, "increasing")
#' 
#' # Generate all 2x2 matrices with decreasing monotonicity
#' matrices_dec <- monotonic_matrices(2, 2, "decreasing")
#'
#' @export
monotonic_matrices <- function(I, J, direction = c("increasing", "decreasing")) {
  direction <- match.arg(direction)
  configs <- grid_upsets(DoseGrid(levels = c(I, J)))
  if (direction == "decreasing") {
    configs <- 1L - configs
  }
  # cells read row by row, 0 before 1
  row_major <- c(t(matrix(seq_len(I * J), I, J)))
  configs <- configs[do.call(order, unname(as.data.frame(configs[, row_major, drop = FALSE]))), ,
                     drop = FALSE]
  lapply(seq_len(nrow(configs)), function(r) matrix(as.numeric(configs[r, ]), I, J))
}
