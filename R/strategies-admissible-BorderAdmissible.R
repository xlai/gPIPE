#' Border Admissibility Rule
#'
#' A reference class that implements an admissibility rule based on the border region
#' between acceptable and unacceptable doses. This rule identifies doses that have
#' a specific number of acceptable and unacceptable neighboring doses.
#'
#' @details
#' The BorderAdmissible class implements the isAdmissible method by identifying dose
#' combinations that lie on the border between acceptable and unacceptable doses.
#' These doses typically have a neighbor sum between 1 and 3, indicating they are
#' at the transition boundary of the dose-toxicity space.
#'
#' @seealso 
#' \code{\link{AdmissibleCombinationRule}} for the base class
#' \code{\link{NeighbourAdmissible}}, \code{\link{ClosestAdmissible}}, 
#' \code{\link{SafetyAdmissible}}, \code{\link{AdjacentAdmissible}} for other admissibility rules
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#'
#' # Create a BorderAdmissible rule
#' border_rule <- BorderAdmissible$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
BorderAdmissible <- setRefClass(
  "BorderAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    #' @description
    #' Determine if dose combinations are admissible based on border position
    #'
    #' @param doseConfig DoseConfiguration object
    #' @param currentDoseLevel Current dose level
    #' @param ... Additional arguments (not used)
    #' @return Vector of indices representing admissible dose combinations
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      border_admissable_level <- findBorderIndices(best_config$currentConfig, n_dose_level)
      return(sort(border_admissable_level))
    }
  )
)

#' Find Border Indices in a Dose Configuration
#'
#' This function identifies indices in a dose configuration that lie on the border
#' between acceptable and unacceptable doses based on neighbor sums.
#'
#' @param A Binary vector representing the current dose configuration
#' @param matrixDims Vector containing the dimensions of the matrix
#'
#' @return Vector of indices representing border positions
#' @keywords internal
findBorderIndices <- function(A, matrixDims) {
  # Step 1: Calculate the neighbour sum
  neighbour_sum <- calculateNeighbourSum(A, matrixDims)
  
  # Step 2: Identify the indices where the neighbour sum is between 1 and 3
  indices <- which(neighbour_sum > 0  & neighbour_sum < 4)

  return(indices)
}
