#' Equal Randomisation Selection Strategy
#'
#' A reference class that implements a selection strategy for choosing the next dose
#' using simple random selection from the admissible doses. This strategy gives equal
#' probability to all admissible doses.
#'
#' @details
#' The EqualRandomisationStrategy class implements the selectDose method by randomly
#' selecting one dose from the set of admissible doses with equal probability. This
#' approach ensures full exploration of the dose space without any bias towards
#' specific doses based on sample sizes or posterior probabilities.
#'
#' @seealso 
#' \code{\link{SelectionStrategy}} for the base class
#' \code{\link{SmallestSampleSizeStrategy}}, \code{\link{PosteriorProbabilityStrategy}},
#' \code{\link{InverseSqrtSampleSizeStrategy}}, \code{\link{InverseDistanceStrategy}} for other selection strategies
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#' 
#' # Create an EqualRandomisationStrategy
#' equal_rand_strategy <- EqualRandomisationStrategy$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
EqualRandomisationStrategy <- setRefClass(
  "EqualRandomisationStrategy",
  contains = "SelectionStrategy",
  methods = list(
    #' @description
    #' Select the next dose using equal randomisation
    #'
    #' @param admissibleDoses Vector of admissible doses
    #' @param summaryStats Summary statistics for each dose combination
    #' @param drugCombi DrugCombinationModel object containing dose information
    #' @return The selected dose combination
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
        selected_level <- sample_one(admissibleDoses)
      return(selected_level)
    }
  )
)