SmallestSampleSizeStrategy <- setRefClass(
  "SmallestSampleSizeStrategy",
  contains = "SelectionStrategy",
  methods = list(
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      sample_size_posterior <- sapply(summaryStats, function(x) x$numPatients) + 
        drugCombi$priorA + drugCombi$priorB
      sample_size_subset <- sample_size_posterior[admissibleDoses]
      selected_level <- sample_one(admissibleDoses[which(sample_size_subset == min(sample_size_subset))])
      return(selected_level)
    }
  )
)