PosteriorProbabilityStrategy <- setRefClass(
  "PosteriorProbabilityStrategy",
  contains = "SelectionStrategy",
  methods = list(
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      probabilities <- drugCombi$p_posterior[admissibleDoses]
      probabilities <- probabilities / sum(probabilities)
      selected_level <- sample_one(admissibleDoses, prob = probabilities)
      return(selected_level)
    }
  )
)