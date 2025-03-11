NeighbourAdmissible <- setRefClass(
  "NeighbourAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      neighbour_level <- findNearestNeighbours(currentDoseLevel, n_dose_level)
      return(neighbour_level)
    }
  )
)


findNearestNeighbours <- function(pointIndex, matrixDims, degree = 1) {
  # Calculate the number of rows and columns from matrixDims
    cols <- matrixDims[1]
    rows <- matrixDims[2]

    # Convert the linear index to row and column indices
    pointRow <- ((pointIndex - 1) %/% cols) + 1
    pointCol <- ((pointIndex - 1) %% cols) + 1
    
    # Create vectors of row and column offsets based on the degree
    rowOffsets <- seq(-degree, degree)
    colOffsets <- seq(-degree, degree)
    
    # Generate a grid of neighbour offsets
    neighbours <- expand.grid(rowOffsets, colOffsets)
    
    # Calculate the row and column indices of neighbours
    neighbourRows <- pointRow + neighbours[,1]
    neighbourCols <- pointCol + neighbours[,2]
    
    # Filter out neighbours that fall outside the matrix boundaries
    validNeighbours <- neighbourRows >= 1 & neighbourRows <= rows &
                        neighbourCols >= 1 & neighbourCols <= cols
    
    # Calculate the linear indices of valid neighbours
    neighbourIndices <- (neighbourRows[validNeighbours] - 1) * cols +
                        neighbourCols[validNeighbours]
    
    # Return the linear indices of the nearest neighbours
    return(neighbourIndices)
}