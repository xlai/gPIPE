#' Inverse Distance Selection Strategy
#'
#' A reference class that implements a selection strategy for choosing the next dose
#' based on the inverse of the distance between posterior probabilities and the target
#' toxicity threshold. This strategy favors doses whose posterior probabilities are
#' closest to the target.
#'
#' @details
#' The InverseDistanceStrategy class implements the selectDose method by calculating
#' probabilities inversely proportional to the distance between each dose's posterior
#' probability and the target toxicity threshold (theta). The distance can be calculated
#' using absolute difference or squared difference.
#'
#' @seealso 
#' \code{\link{SelectionStrategy}} for the base class
#' \code{\link{SmallestSampleSizeStrategy}}, \code{\link{PosteriorProbabilityStrategy}},
#' \code{\link{InverseSqrtSampleSizeStrategy}}, \code{\link{EqualRandomisationStrategy}} for other selection strategies
#'
#' @examples
#' \dontrun{
#' # Create base admissible rule class first if needed
#' library(gPIPE)
#' 
#' # Create an InverseDistanceStrategy using absolute distance
#' inverse_dist_strategy <- InverseDistanceStrategy$new()
#' }
#'
#' @importFrom methods setRefClass new
#' @export
InverseDistanceStrategy <- setRefClass(
  "InverseDistanceStrategy",
  contains = "SelectionStrategy",
  methods = list(
    #' @description
    #' Select the next dose based on inverse distance to target threshold
    #'
    #' @param admissibleDoses Vector of admissible doses
    #' @param summaryStats Summary statistics for each dose combination
    #' @param drugCombi DrugCombinationModel object containing posterior probabilities and target threshold
    #' @param distance_type Type of distance metric to use: "absolute" or "squared"
    #' @return The selected dose combination
    selectDose = function(admissibleDoses, summaryStats, drugCombi, distance_type = "absolute") {
      p_posterior <- drugCombi$p_posterior[admissibleDoses]
      dist_tolerance <- 0.01
      distances <- if (distance_type == "absolute") {
        abs(p_posterior - drugCombi$theta + dist_tolerance)
      } else if (distance_type == "squared") {
        (p_posterior - drugCombi$theta + dist_tolerance) ^ 2
      } else {
        stop("Invalid distance_type. Use 'absolute' or 'squared'.")
      }
      probabilities <- 1 / distances
      probabilities <- probabilities / sum(probabilities)
      selected_level <- sample_one(admissibleDoses, prob = probabilities)
      return(selected_level)
    }
  )
)
