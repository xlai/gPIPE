DrugCombi <- setRefClass(
    "DrugCombination",
    fields = list(
        drugs = "list",  # List to store Drug objects
        doseCombinations = "list"  # List to store generated combinations
    ),
    methods = list(
        initialize = function(drugs = list()) {
            drugs <<- drugs
            doseCombinations <<- list()
            .self$generateCombinations()
        },
        addDrug = function(drug) {
            if (!inherits(drug, "Drug")) {
                stop("Invalid object: The drug must be a Drug object.")
            }
            drugs <<- c(drugs, list(drug))
            .self$generateCombinations()
        },
        generateCombinations = function() {
            # Generate all valid combinations of dose levels
            if (length(drugs) > 0) {
                combinationLists <- lapply(drugs, function(d) d$getDoseLevels())
                doseCombinations <<- expand.grid(combinationLists)
                # Validate and filter combinations here if needed
            }
        },
        getDoseCombinations = function() {
            return(doseCombinations)
        },
        print = function() {
            cat("Drug Combination:\n")
            print(doseCombinations)
        }
        # ... (other methods as needed)
    )
)
