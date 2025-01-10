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

findBorderIndices <- function(A, matrixDims) {
  # Step 1: Calculate the neighbour sum
  neighbour_sum <- calculateNeighbourSum(A, matrixDims)
  
  # Step 2: Identify the indices where the neighbour sum is between 1 and 3
  indices <- which(neighbour_sum > 1  & neighbour_sum < 3)

  return(indices)
}
