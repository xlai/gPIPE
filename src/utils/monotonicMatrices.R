monotonic_matrices <- function(I, J, direction = c("increasing", "decreasing")) {
  direction <- match.arg(direction)
  
  # Function to check if a matrix is monotonic
  is_monotonic <- function(m, dir) {
    # Check rows
    for (r in 1:nrow(m)) {
      row_diffs <- diff(m[r,])
      if ((dir == "decreasing" && any(row_diffs > 0)) || 
          (dir == "increasing" && any(row_diffs < 0))) {
        return(FALSE)
      }
    }
    
    # Check columns
    for (c in 1:ncol(m)) {
      col_diffs <- diff(m[,c])
      if ((dir == "decreasing" && any(col_diffs > 0)) || 
          (dir == "increasing" && any(col_diffs < 0))) {
        return(FALSE)
      }
    }
    
    return(TRUE)
  }
  
  # Recursive function to build matrices
  build_matrix <- function(m, i, j) {
    # If matrix is complete, return it
    if (i > I) {
      return(list(m))
    }
    
    # Move to next position
    next_i <- i
    next_j <- j + 1
    if (next_j > J) {
      next_i <- i + 1
      next_j <- 1
    }
    
    # Determine valid values for current position
    valid_values <- c(0, 1)
    
    # Apply constraints based on neighbors
    if (direction == "increasing") {
      if (j > 1) {
        valid_values <- valid_values[valid_values >= m[i, j-1]]
      }
      if (i > 1) {
        valid_values <- valid_values[valid_values >= m[i-1, j]]
      }
    } else {
      if (j > 1) {
        valid_values <- valid_values[valid_values <= m[i, j-1]]
      }
      if (i > 1) {
        valid_values <- valid_values[valid_values <= m[i-1, j]]
      }
    }
    
    # Try each valid value
    result <- list()
    for (val in valid_values) {
      m[i, j] <- val
      result <- c(result, build_matrix(m, next_i, next_j))
    }
    
    return(result)
  }
  
  # Start with an empty matrix
  empty_matrix <- matrix(NA, nrow = I, ncol = J)
  return(build_matrix(empty_matrix, 1, 1))
}
