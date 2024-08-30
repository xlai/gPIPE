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
            labelCombinations <- labelCombinations[order(labelCombinations$Var1, labelCombinations$Var2),] # reorder
            # Creating combination names and mapping to numeric levels
            combinedNames <- apply(labelCombinations, 1, paste, collapse = ".")

            # Mapping each label combination to its corresponding numeric levels
            levelCombinations <- expand.grid(levelsList, stringsAsFactors = FALSE)
            levelCombinations <- levelCombinations[order(levelCombinations$Var1, levelCombinations$Var2),] # reorder

            # Assigning row numbers as unique identifiers for each combination
            combinationIDs <- seq_len(nrow(levelCombinations))
            
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
        getNumberOfDoseLevels = function(combined = FALSE) {
            if (combined) {
                # Return the number of combined dose levels
                return(length(doseCombinations))
            } else {
                # Return a tuple of dose levels for each separate drug
                return(sapply(drugs, function(drug) drug$getNumberOfDoseLevels()))
            }
        },
        getDoseCombinationsLevel = function(doseCombi_selected = NULL) {
            if (is.null(doseCombi_selected)){
                return(sapply(doseCombinations, function(x) x$id))
            }
            else {
               return(sapply(doseCombinations[doseCombi_selected], function(x) x$id))
            }
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
