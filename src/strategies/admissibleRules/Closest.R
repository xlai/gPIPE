ClosestAdmissible <- setRefClass(
  "ClosestAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      closest_admissable_level <- calculateClosestDoses(best_config$currentConfig, n_dose_level[1], n_dose_level[2])
      return(which(closest_admissable_level))
    }
  )
)

calculateClosestDoses <- function(vec, nrows, ncols) {
  # Convert the vector into a matrix
  adjacent_doses <- calculateAdjacentDoses(vec, nrows, ncols)
  vec_matrix <- matrix(vec, nrows, ncols, byrow = FALSE)
  # Identifying dominant upper and lower doses
	dominantu<- adjacent_doses * vec_matrix * (rbind(0,vec_matrix[-nrows,]) == 0 & cbind(0,vec_matrix[,-ncols]) == 0)
	dominantl<- adjacent_doses * (1 - vec_matrix) * (rbind(vec_matrix[-1,], 1) == 1 & cbind(vec_matrix[,-1],1) == 1)
	dominant<-dominantl | dominantu
  
  return(dominant)
}