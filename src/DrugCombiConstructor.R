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
            if (length(drugs) == 0) {
                doseCombinations <<- list()
                return()
            }

            # Extracting dose labels and levels from each drug
            labelsList <- lapply(drugs, function(d) d$getDoseLabels())
            levelsList <- lapply(drugs, function(d) d$getDoseLevels())

            # Creating all possible combinations of dose labels
            labelCombinations <- expand.grid(labelsList, stringsAsFactors = FALSE)
            
            # Creating combination names and mapping to numeric levels
            combinedNames <- apply(labelCombinations, 1, paste, collapse = "-")

            # Assigning row numbers as unique identifiers for each combination
            combinationIDs <- seq_len(nrow(levelCombinations))

            # Mapping each label combination to its corresponding numeric levels
            levelCombinations <- expand.grid(levelsList, stringsAsFactors = FALSE)
            combinedLevels <- apply(levelCombinations, 1, function(x) as.numeric(x))

            # Storing combinations in doseCombinations
            # Each combination is a list with levels and its corresponding combination ID
            # Storing combinations in doseCombinations
            doseCombinations <<- setNames(
                lapply(combinationIDs, function(i) {
                    list(levels = combinedLevels[, i], id = i)
                }),
                combinedNames
            )
        },
        updateCombinations = function() {
            # Call this method to update the combinations
            # whenever there's a change in any of the Drug objects
            .self$generateCombinations()
        },        
        getDoseCombinations = function() {
            return(doseCombinations)
        },
        print = function() {
            cat("Drug Combination of:")
            drugNames <- lapply(drugs, function(drug) drug$name)
            cat(paste(unlist(drugNames), collapse = ", "), "\n")

            cat("Drug Levels:\n")
            if (length(doseCombinations) == 0) {
                cat("No dose combinations generated yet.\n")
            } else {
                    cat( names(doseCombinations), "\n")
            }
        }
        # ... (other methods as needed)
    )
)
