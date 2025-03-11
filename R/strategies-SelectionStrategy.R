#' Dose Selection Strategy Base Class
#'
#' An abstract reference class that defines the interface for selecting the next
#' dose combination in a dose-finding study. This class serves as a base class
#' for implementing specific dose selection strategies.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{initialize(...)}}{Initialize a new SelectionStrategy object}
#'   \item{\code{selectDose(admissibleDoses, summaryStats, drugCombi)}}{Select the next dose from admissible doses}
#' }
#'
#' @details
#' The SelectionStrategy class provides a framework for implementing different strategies
#' to select the next dose combination to be administered in a dose-finding study.
#' This selection is typically made from the set of admissible doses determined by
#' an AdmissibleCombinationRule object.
#'
#' @seealso \code{\link{AdmissibleCombinationRule}} for the complementary rule class that determines
#' which doses are admissible
#' @importFrom methods setRefClass new
#' @export
SelectionStrategy <- setRefClass(
  "SelectionStrategy",
  methods = list(
    #' @description
    #' Initialize a new SelectionStrategy object
    initialize = function(...) {
      cat("Initializing SelectionStrategy object...\n")
    },    
    
    #' @description
    #' Select the next dose from admissible doses
    #'
    #' @param admissibleDoses Vector of admissible doses
    #' @param summaryStats Summary statistics for each dose combination
    #' @param drugCombi DrugCombination object containing dose information
    #' @return The selected dose combination
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      # Placeholder: Implement in subclasses
      stop("Method 'selectDose' must be implemented by subclasses")
    }
  )
)

#' Sample a Single Element from a Vector
#'
#' This helper function safely samples a single element from a vector,
#' handling edge cases such as empty vectors or vectors with only one element.
#'
#' @param vec Vector to sample from
#' @param prob Optional vector of probability weights for sampling
#'
#' @return A single sampled element from the vector, or NA if the vector is empty
#' @keywords internal
sample_one <- function(vec, prob = NULL) {
    if (length(vec) == 0) {
        return(NA)  # Or any appropriate default value
    } else if (length(vec) == 1) {
        return(vec)  # Return the single element directly
    } else {
        return(sample(vec, size = 1, prob = prob))
    }
}
