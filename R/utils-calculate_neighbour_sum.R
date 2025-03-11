#' Calculate Sum of Neighboring Elements in a Matrix
#'
#' This function calculates the sum of adjacent elements (up, down, left, and right)
#' for each element in a matrix. The function is particularly useful in dose-finding studies
#' for identifying boundary regions between acceptable and unacceptable dose combinations.
#'
#' @param A A numeric vector representing a flattened matrix
#' @param matrixDims A numeric vector of length 2 containing the dimensions of the matrix:
#'        c(number_of_columns, number_of_rows)
#'
#' @return A matrix of the same dimensions as the input matrix, where each element
#'         contains the sum of its neighboring elements (up, down, left, and right)
#'
#' @details
#' The function reshapes the input vector into a matrix of dimensions specified by `matrixDims`,
#' pads the matrix with a border (1s on top, left, and part of right/bottom, 0s on the 
#' bottom-right corner), and calculates the sum of the four adjacent elements for each position.
#' 
#' This is particularly useful for identifying transition boundaries in dose-toxicity matrices
#' where a sum of 2 often indicates potential Maximum Tolerated Dose (MTD) positions.
#'
#' @examples
#' # Create a sample binary matrix
#' mat <- c(1, 1, 0, 1, 0, 0, 0, 0, 0)
#' dims <- c(3, 3)  # 3 columns, 3 rows
#' 
#' # Calculate neighbor sums
#' neighbor_sums <- calculateNeighbourSum(mat, dims)
#' 
#' @export
calculateNeighbourSum <- function(A, matrixDims) {

  # Calculate the number of rows and columns from matrixDims
    ncols <- matrixDims[1]
    nrows <- matrixDims[2] 
  
  # Reshape the vector back to a matrix
  mat <- matrix(A, nrow = nrows, ncol = ncols)
  
  # Create a padded matrix with extra row and column at the borders
  padded_A <- matrix(1, nrow = nrows + 2, ncol = ncols + 2)
  padded_A[nrows + 2,] <- 0
  padded_A[,ncols + 2] <- 0

  # Place the original matrix in the centre of the padded matrix
  padded_A[2:(nrows + 1), 2:(ncols + 1)] <- mat
  
  # Calculate the sum of the neighbours (up, down, left, right)
  neighbour_sum <- padded_A[1:nrows, 2:(ncols + 1)] +  # up
                   padded_A[3:(nrows + 2), 2:(ncols + 1)] +  # down
                   padded_A[2:(nrows + 1), 1:ncols] +  # left
                   padded_A[2:(nrows + 1), 3:(ncols + 2)]    # right
  
  return(neighbour_sum)
}
