InverseDistanceStrategy <- setRefClass(
  "InverseDistanceStrategy",
  contains = "SelectionStrategy",
  methods = list(
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
