#' Adjacent Admissibility Rule
#'
#' A reference class that implements an admissibility rule based on adjacency to the
#' current set of acceptable and unacceptable doses. This rule identifies doses that
#' are immediately adjacent to the acceptance boundary.
#'
#' @details
#' The AdjacentAdmissible class implements the isAdmissible method by identifying dose
#' combinations that are adjacent to the current boundary between acceptable and unacceptable
#' doses. This approach helps in efficiently exploring the dose-toxicity relationship by
#' focusing on the region where the dose-toxicity relationship is changing.
#'
#' @seealso 
#' \code{\link{AdmissibleCombinationRule}} for the base class
#' \code{\link{NeighbourAdmissible}}, \code{\link{ClosestAdmissible}}, 
#' \code{\link{BorderAdmissible}}, \code{\link{SafetyAdmissible}} for other admissibility rules
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#' 
#' # Create an AdjacentAdmissible rule
#' adjacent_rule <- AdjacentAdmissible$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
AdjacentAdmissible <- setRefClass(
  "AdjacentAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    #' @description
    #' Determine if dose combinations are admissible based on adjacency to the acceptance boundary
    #'
    #' @param doseConfig DoseConfiguration object
    #' @param currentDoseLevel Current dose level
    #' @param ... Additional arguments (not used)
    #' @return Vector of indices representing admissible dose combinations
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      adjacent_admissable_level <- calculateAdjacentDoses(doseConfig$currentConfig, n_dose_level[1], n_dose_level[2])
      return(c(as.numeric(adjacent_admissable_level)))
    }
  )
)

#' Calculate Adjacent Doses
#'
#' This function identifies dose combinations that are adjacent to the boundary
#' between acceptable and unacceptable doses.
#'
#' @param vec Binary vector representing the current dose configuration
#' @param nrows Number of rows in the dose matrix
#' @param ncols Number of columns in the dose matrix
#'
#' @return Binary matrix where TRUE indicates adjacent doses
#' @keywords internal
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