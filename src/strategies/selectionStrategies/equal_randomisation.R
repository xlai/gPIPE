EqualRandomisationStrategy <- setRefClass(
  "EqualRandomisationStrategy",
  contains = "SelectionStrategy",
  methods = list(
    selectDose = function(admissibleDoses, summaryStats, drugCombi) {
        selected_level <- sample_one(admissibleDoses)
      return(selected_level)
    }
  )
)