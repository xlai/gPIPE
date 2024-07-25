SelectionStrategy <- setRefClass(
  "SelectionStrategy",
  methods = list(
    initialize = function(...) {
      cat("Initializing SelectionStrategy object...\n")
    },    
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      # implement selecting one(s) with smallest sample size
      sample_size_posterior <- sapply(summaryStats, function(x) x$numPatients) + 
        drugCombi$priorA + drugCombi$priorB
      # Find indices in 'b' that correspond to the smallest value in 'a' and randomly select one
      sample_size_subset <- sample_size_posterior[admissibleDoses]
      selected_level <- sample(admissibleDoses[which(sample_size_subset == min(sample_size_subset))], 1)
      return(selected_level)
    }
  )
)