BorderAdmissible <- setRefClass(
  "BorderAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      border_admissable_level <- findBorderIndices(best_config$currentConfig, n_dose_level)
      return(sort(border_admissable_level))
    }
  )
)

findBorderIndices <- function(vec_mat, matrixDims) {
  # Calculate the number of rows and columns from matrixDims
    ncols <- matrixDims[1]
    nrows <- matrixDims[2]  
  # Reshape the vector back to a matrix
  mat <- matrix(vec_mat, nrow = nrows, ncol = ncols, byrow = TRUE)
  
  # Function to convert (i, j) to linear index
  linear_index <- function(i, j, ncols) {
    return((i - 1) * ncols + j)
  }
  
  # Identify row-wise transitions from 0 to 1
  row_transitions <- which(mat[, -ncols] == 0 & mat[, -1] == 1, arr.ind = TRUE)
  row_indices <- unique(c(
    linear_index(row_transitions[, 1], row_transitions[, 2], ncols),
    linear_index(row_transitions[, 1], row_transitions[, 2] + 1, ncols)
  ))
  
  # Identify column-wise transitions from 0 to 1
  col_transitions <- which(mat[-nrows, ] == 0 & mat[-1, ] == 1, arr.ind = TRUE)
  col_indices <- unique(c(
    linear_index(col_transitions[, 1], col_transitions[, 2], ncols),
    linear_index(col_transitions[, 1] + 1, col_transitions[, 2], ncols)
  ))
  
  # Combine and return unique indices
  border_indices <- unique(c(row_indices, col_indices))
  
  return(border_indices)
}