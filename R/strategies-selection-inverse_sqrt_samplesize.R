#' Inverse Square Root Sample Size Selection Strategy
#'
#' A reference class that implements a selection strategy for choosing the next dose
#' based on the inverse of the square root of sample sizes. This strategy balances
#' exploration of less-tested doses with a more moderate weighting than the simple
#' inverse sample size approach.
#'
#' @details
#' The InverseSqrtSampleSizeStrategy class implements the selectDose method by calculating
#' probabilities inversely proportional to the square root of the posterior sample size
#' (number of patients plus prior pseudocounts) for each admissible dose. This approach
#' favors doses with fewer observations but with a less extreme bias than using raw
#' inverse sample sizes.
#'
#' @seealso 
#' \code{\link{SelectionStrategy}} for the base class
#' \code{\link{SmallestSampleSizeStrategy}}, \code{\link{PosteriorProbabilityStrategy}},
#' \code{\link{InverseDistanceStrategy}}, \code{\link{EqualRandomisationStrategy}} for other selection strategies
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#' 
#' # Create an InverseSqrtSampleSizeStrategy
#' inverse_sqrt_strategy <- InverseSqrtSampleSizeStrategy$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
InverseSqrtSampleSizeStrategy <- setRefClass(
  "InverseSqrtSampleSizeStrategy",
  contains = "SelectionStrategy",
  methods = list(
    #' @description
    #' Select the next dose based on the inverse square root of sample sizes
    #'
    #' @param admissibleDoses Vector of admissible doses
    #' @param summaryStats Summary statistics for each dose combination
    #' @param drugCombi DrugCombinationModel object containing dose information
    #' @return The selected dose combination
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      # Calculate the posterior sample size for each dose
      sample_size_posterior <- sapply(summaryStats, function(x) x$numPatients) + 
        drugCombi$priorA + drugCombi$priorB
      
      # Subset to only include admissible doses
      sample_size_subset <- sample_size_posterior[admissibleDoses]
      
      # Calculate probabilities inversely proportional to the square root of sample sizes
      probabilities <- 1 / sqrt(sample_size_subset)
      
      # Normalize the probabilities
      probabilities <- probabilities / sum(probabilities)
      
      # Select a dose based on the normalized probabilities
      selected_level <- sample_one(admissibleDoses, prob = probabilities)
      
      return(selected_level)
    }
  )
)
