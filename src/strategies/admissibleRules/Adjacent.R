AdjacentAdmissible <- setRefClass(
  "AdjacentAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      adjacent_admissable_level <- calculateAdjacentDoses(doseConfig$currentConfig, n_dose_level[1], n_dose_level[2])
      return(c(as.numeric(adjacent_admissable_level)))
    }
  )
)

calculateAdjacentDoses <- function(vec, nrows, ncols) {
  # Ensure there are more than one level for both drugs
  if(nrows < 2 || ncols < 2) {
    stop("Admissible doses can only be calculated when both drugs have more than one level")
  }
  
  # Convert the vector into a matrix
  mat <- matrix(vec, nrow = nrows, ncol = ncols, byrow = FALSE)
  
  # Identifying dominant upper and lower doses
  dominantu <- mat == 1 & (rbind(0, mat[-nrows,]) == 0 | cbind(0, mat[,-ncols]) == 0 | rbind(0, cbind(0, mat[,-ncols])[-nrows,]) == 0)
  dominantl <- mat == 0 & (rbind(mat[-1,], 1) == 1 | cbind(mat[,-1], 1) == 1 | rbind(cbind(mat[,-1], 1)[-1,], 1) == 1)
  dominant <- dominantl | dominantu
  
  return(dominant)
}