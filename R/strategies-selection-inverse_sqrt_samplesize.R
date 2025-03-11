InverseSqrtSampleSizeStrategy <- setRefClass(
  "InverseSqrtSampleSizeStrategy",
  contains = "SelectionStrategy",
  methods = list(
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
