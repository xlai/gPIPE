SelectionStrategy <- setRefClass(
  "SelectionStrategy",
  methods = list(
    initialize = function(...) {
      cat("Initializing SelectionStrategy object...\n")
    },    
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
      # Placeholder: Implement in subclasses
    }
  )
)

sample_one <- function(vec, prob = NULL) {
    if (length(vec) == 0) {
        return(NA)  # Or any appropriate default value
    } else if (length(vec) == 1) {
        return(vec)  # Return the single element directly
    } else {
        return(sample(vec, size = 1, prob = prob))
    }
}
