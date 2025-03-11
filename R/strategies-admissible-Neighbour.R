#' Neighbour-Based Admissibility Rule
#'
#' A reference class that implements a neighbour-based rule for determining which dose combinations
#' are admissible in a dose-finding study. This rule restricts dose escalation or de-escalation
#' to neighbouring dose combinations.
#'
#' @details
#' The NeighbourAdmissible class implements the isAdmissible method by identifying dose combinations
#' that are adjacent to the current dose level. This ensures that dose changes occur gradually,
#' enhancing patient safety by preventing large jumps in dosing.
#'
#' @seealso 
#' \code{\link{AdmissibleCombinationRule}} for the base class
#' \code{\link{SafetyAdmissible}}, \code{\link{ClosestAdmissible}}, 
#' \code{\link{BorderAdmissible}}, \code{\link{AdjacentAdmissible}} for other admissibility rules
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#'
#' # Create a NeighbourAdmissible rule
#' neighbour_rule <- NeighbourAdmissible$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
NeighbourAdmissible <- setRefClass(
  "NeighbourAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    #' @description
    #' Determine if dose combinations are admissible based on proximity to current dose
    #'
    #' @param doseConfig DoseConfiguration object
    #' @param currentDoseLevel Current dose level
    #' @param ... Additional arguments (not used)
    #' @return Vector of indices representing admissible dose combinations
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      neighbour_level <- findNearestNeighbours(currentDoseLevel, n_dose_level)
      return(neighbour_level)
    }
  )
)

#' Find Nearest Neighbours in a Matrix
#'
#' This function identifies the indices of neighboring cells in a matrix given a point index.
#'
#' @param pointIndex Integer index of the reference point in the matrix
#' @param matrixDims Vector containing the dimensions of the matrix (columns, rows)
#' @param degree Integer indicating the neighborhood degree (1 for immediate neighbors)
#'
#' @return Vector of indices representing neighbors of the reference point
#' @keywords internal
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