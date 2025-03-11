#' Closest Boundary Admissibility Rule
#'
#' A reference class that implements an admissibility rule based on proximity to
#' the boundary between acceptable and unacceptable doses. This rule identifies
#' doses that are closest to the efficacy-toxicity threshold.
#'
#' @details
#' The ClosestAdmissible class implements the isAdmissible method by identifying dose
#' combinations that lie on the boundary between acceptable and unacceptable doses.
#' These doses are closest to the maximum tolerated dose (MTD) and represent the most
#' informative doses to test next.
#'
#' @seealso 
#' \code{\link{AdmissibleCombinationRule}} for the base class
#' \code{\link{NeighbourAdmissible}}, \code{\link{SafetyAdmissible}}, 
#' \code{\link{BorderAdmissible}}, \code{\link{AdjacentAdmissible}} for other admissibility rules
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#'
#' # Create a ClosestAdmissible rule
#' closest_rule <- ClosestAdmissible$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
ClosestAdmissible <- setRefClass(
  "ClosestAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    #' @description
    #' Determine if dose combinations are admissible based on proximity to the efficacy-toxicity boundary
    #'
    #' @param doseConfig DoseConfiguration object
    #' @param currentDoseLevel Current dose level
    #' @param ... Additional arguments (not used)
    #' @return Vector of indices representing admissible dose combinations
    isAdmissible = function(doseConfig, currentDoseLevel, ...) {
      # Custom implementation
      best_config <- doseConfig$bestConfigs
      n_dose_level <- best_config$drugCombi$getNumberOfDoseLevels()
      closest_admissable_level <- calculateClosestDoses(best_config$currentConfig, n_dose_level[1], n_dose_level[2])
      return(which(closest_admissable_level))
    }
  )
)

#' Calculate Doses Closest to the Boundary
#'
#' This function identifies dose combinations that lie on the boundary between
#' acceptable and unacceptable doses.
#'
#' @param vec Binary vector representing the current dose configuration
#' @param nrows Number of rows in the dose matrix
#' @param ncols Number of columns in the dose matrix
#'
#' @return Binary matrix where TRUE indicates doses closest to the boundary
#' @keywords internal
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