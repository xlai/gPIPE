#' Safety-Based Admissibility Rule
#'
#' A reference class that implements a safety-based rule for determining which dose combinations
#' are admissible in a dose-finding study. This rule evaluates dose combinations based on
#' a safety threshold applied to the posterior probabilities.
#'
#' @details
#' The SafetyAdmissible class implements the isAdmissible method by comparing the normalized
#' posterior probabilities against a safety threshold. Doses with posterior probabilities 
#' below the threshold are considered admissible.
#'
#' @seealso 
#' \code{\link{AdmissibleCombinationRule}} for the base class
#' \code{\link{NeighbourAdmissible}}, \code{\link{ClosestAdmissible}}, 
#' \code{\link{BorderAdmissible}}, \code{\link{AdjacentAdmissible}} for other admissibility rules
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#'
#' # Create a SafetyAdmissible rule with default safety threshold
#' safety_rule <- SafetyAdmissible$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
SafetyAdmissible <- setRefClass(
  "SafetyAdmissible",
  contains = "AdmissibleCombinationRule", 
  methods = list(
    #' @description
    #' Determine if dose combinations are admissible based on safety threshold
    #'
    #' @param doseConfig DoseConfiguration object
    #' @param currentDoseLevel Current dose level
    #' @param drugCombiModel Drug combination model with posterior probabilities
    #' @param safety_threshold Threshold for safety (default: 0.8)
    #' @return Vector of indices of admissible dose combinations
    isAdmissible = function(doseConfig, currentDoseLevel, drugCombiModel, safety_threshold = 0.8) {
      # Custom implementation
      doseConfig_posterior <- doseConfig$updatePipeEstimator(drugCombiModel$p_posterior)
      sum_posterior <- Reduce("+", doseConfig_posterior)
      normalized_posterior <- lapply(doseConfig_posterior, function(x) x / sum_posterior)

      validConfig_list <- lapply(doseConfig$validConfigs, function(x) x$currentConfig)

      # Use Map to multiply element-wise
      result <- Map(function(x, y) x * y, normalized_posterior, validConfig_list)
      result_collapsed <-  Reduce("+", result)

      return(which(result_collapsed < safety_threshold))
    }
  )
)