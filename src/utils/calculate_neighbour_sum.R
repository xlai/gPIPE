calculateNeighbourSum <- function(A, matrixDims) {

  # This function takes a vectorized matrix A and its dimensions (matrixDims) as input,
  # reshapes A back into its matrix form, pads the matrix with a border, and then calculates
  # the sum of the neighbouring elements (up, down, left, and right) for each element in the matrix.
  # The result is a matrix of the same size where each element contains the sum of its neighbours.
      
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
