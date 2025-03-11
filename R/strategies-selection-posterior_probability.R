#' Posterior Probability Selection Strategy
#'
#' A reference class that implements a selection strategy for choosing the next dose
#' based on posterior probabilities. This strategy selects doses with probabilities
#' proportional to their posterior probabilities of being the target dose.
#'
#' @details
#' The PosteriorProbabilityStrategy class implements the selectDose method by using
#' the posterior probabilities from the drug combination model to weight the selection
#' of admissible doses. Doses with higher posterior probabilities have a greater chance
#' of being selected.
#'
#' @seealso 
#' \code{\link{SelectionStrategy}} for the base class
#' \code{\link{SmallestSampleSizeStrategy}}, \code{\link{InverseSqrtSampleSizeStrategy}},
#' \code{\link{InverseDistanceStrategy}}, \code{\link{EqualRandomisationStrategy}} for other selection strategies
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#'
#' # Create a PosteriorProbabilityStrategy
#' posterior_prob_strategy <- PosteriorProbabilityStrategy$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
PosteriorProbabilityStrategy <- setRefClass(
  "PosteriorProbabilityStrategy",
  contains = "SelectionStrategy",
  methods = list(
    #' @description
    #' Select the next dose based on posterior probabilities
    #'
    #' @param admissibleDoses Vector of admissible doses
    #' @param summaryStats Summary statistics for each dose combination
    #' @param drugCombi DrugCombinationModel object containing posterior probabilities
    #' @return The selected dose combination
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      probabilities <- drugCombi$p_posterior[admissibleDoses]
      probabilities <- probabilities / sum(probabilities)
      selected_level <- sample_one(admissibleDoses, prob = probabilities)
      return(selected_level)
    }
  )
)