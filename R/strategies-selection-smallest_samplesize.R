#' Smallest Sample Size Selection Strategy
#'
#' A reference class that implements a selection strategy for choosing the next dose
#' based on the smallest sample size. This strategy prioritizes dose levels that have
#' been tested on the fewest patients, promoting balanced exploration of the dose space.
#'
#' @details
#' The SmallestSampleSizeStrategy class implements the selectDose method by identifying
#' the admissible doses with the smallest posterior sample size (number of patients plus
#' prior pseudocounts). If multiple doses have the same minimum sample size, one is 
#' randomly selected.
#'
#' @seealso 
#' \code{\link{SelectionStrategy}} for the base class
#' \code{\link{PosteriorProbabilityStrategy}}, \code{\link{InverseSqrtSampleSizeStrategy}},
#' \code{\link{InverseDistanceStrategy}}, \code{\link{EqualRandomisationStrategy}} for other selection strategies
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#'
#' # Create a SmallestSampleSizeStrategy
#' smallest_sample_strategy <- SmallestSampleSizeStrategy$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
SmallestSampleSizeStrategy <- setRefClass(
  "SmallestSampleSizeStrategy",
  contains = "SelectionStrategy",
  methods = list(
    #' @description
    #' Select the next dose based on the smallest sample size
    #'
    #' @param admissibleDoses Vector of admissible doses
    #' @param summaryStats Summary statistics for each dose combination
    #' @param drugCombi DrugCombinationModel object containing dose information
    #' @return The selected dose combination
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      sample_size_posterior <- sapply(summaryStats, function(x) x$numPatients) + 
        drugCombi$priorA + drugCombi$priorB
      sample_size_subset <- sample_size_posterior[admissibleDoses]
      selected_level <- sample_one(admissibleDoses[which(sample_size_subset == min(sample_size_subset))])
      return(selected_level)
    }
  )
)