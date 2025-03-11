#' Drug Combination Class
#'
#' A reference class that represents combinations of drugs.
#' This class manages objects of the \code{\link{Drug-class}} class to create
#' and manage combinations of different drugs at various dose levels.
#'
#' @field drugs List of Drug objects
#' @field doseCombinations List of generated dose combinations
#'
#' @seealso \code{\link{Drug-class}} for the base drug class
#' @importFrom methods setRefClass new
#' @export
DrugCombi <- setRefClass(
    "DrugCombination",
    fields = list(
        drugs = "list",  # List to store Drug objects
        doseCombinations = "list"  # List to store generated combinations
    ),
    methods = list(
        #' @description
        #' Initialize a new DrugCombination object
        #'
        #' @param drugs List of Drug objects
        initialize = function(drugs = list()) {
            drugs <<- drugs
            doseCombinations <<- list()
            .self$generateCombinations()
        },
        
        #' @description
        #' Add a new drug to the combination
        #'
        #' @param drug A Drug object to add
        addDrug = function(drug) {
            if (!inherits(drug, "Drug")) {
                stop("Invalid object: The drug must be a Drug object.")
            }
            drugs <<- c(drugs, list(drug))
            .self$generateCombinations()
        },
        
        #' @description
        #' Generate all possible combinations of drug dose levels
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
            doseCombinations <<- setNames(
                lapply(combinationIDs, function(i) {
                    list(levels = combinedLevels[, i], id = i)
                }),
                combinedNames
            )
        },
        
        #' @description
        #' Update the combinations whenever there's a change in any of the Drug objects
        updateCombinations = function() {
            .self$generateCombinations()
        },
        
        #' @description
        #' Get all dose combinations
        #'
        #' @return List of dose combinations
        getDoseCombinations = function() {
            return(doseCombinations)
        },
        
        #' @description
        #' Get the number of dose levels
        #'
        #' @param combined Logical indicating whether to return the total number of combined dose levels (TRUE) or a vector of dose levels for each drug (FALSE)
        #' @return Numeric value or vector of the number of dose levels
        getNumberOfDoseLevels = function(combined = FALSE) {
            if (combined) {
                # Return the number of combined dose levels
                return(length(doseCombinations))
            } else {
                # Return a tuple of dose levels for each separate drug
                return(sapply(drugs, function(drug) drug$getNumberOfDoseLevels()))
            }
        },
        
        #' @description
        #' Get the level identifiers for dose combinations
        #'
        #' @param doseCombi_selected Optional parameter to specify which dose combinations to return levels for
        #' @return Vector of dose combination level identifiers
        getDoseCombinationsLevel = function(doseCombi_selected = NULL) {
            if (is.null(doseCombi_selected)){
                return(sapply(doseCombinations, function(x) x$id))
            }
            else {
               return(sapply(doseCombinations[doseCombi_selected], function(x) x$id))
            }
        },
        
        #' @description
        #' Print a summary of the drug combination
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
    )
)

#' Create a Drug Combination object
#'
#' This function creates a new DrugCombi object for managing combinations of drugs
#' in dose-finding studies. It handles simple initialization with parameters or with
#' pre-existing Drug objects.
#'
#' @param drugs List of \code{\link{Drug-class}} objects
#' @param drug_A Optional list with settings for first drug, containing 'name' and 'dose_levels'
#' @param drug_B Optional list with settings for second drug, containing 'name' and 'dose_levels'
#'
#' @return A new DrugCombi object
#' 
#' @details
#' This function provides a convenient way to create a DrugCombi object either from
#' existing Drug objects or by specifying the parameters for two drugs directly.
#' 
#' If `drugs` is provided, it should be a list of Drug objects.
#' If `drug_A` and `drug_B` are provided, new Drug objects will be created with the specified parameters.
#'
#' @seealso \code{\link{Drug-class}} for the drug class, \code{\link{createDrug}} for creating Drug objects
#' @export
createDrugCombi <- function(drugs = NULL, drug_A = NULL, drug_B = NULL) {
  # If the drugs list is provided, use it directly
  if (!is.null(drugs)) {
    return(DrugCombi$new(drugs = drugs))
  }
  
  # If drug_A and drug_B parameters are provided, create Drug objects
  drugs_list <- list()
  
  if (!is.null(drug_A)) {
    if (is.list(drug_A) && "name" %in% names(drug_A)) {
      if ("dose_levels" %in% names(drug_A)) {
        drug1 <- Drug$new(name = drug_A$name, dose_levels = drug_A$dose_levels)
      } else if ("doseCount" %in% names(drug_A)) {
        drug1 <- createDrug(name = drug_A$name, doseCount = drug_A$doseCount)
      } else {
        drug1 <- createDrug(name = drug_A$name, doseCount = 3)
      }
      drugs_list <- c(drugs_list, list(drug1))
    }
  }
  
  if (!is.null(drug_B)) {
    if (is.list(drug_B) && "name" %in% names(drug_B)) {
      if ("dose_levels" %in% names(drug_B)) {
        drug2 <- Drug$new(name = drug_B$name, dose_levels = drug_B$dose_levels)
      } else if ("doseCount" %in% names(drug_B)) {
        drug2 <- createDrug(name = drug_B$name, doseCount = drug_B$doseCount)
      } else {
        drug2 <- createDrug(name = drug_B$name, doseCount = 3)
      }
      drugs_list <- c(drugs_list, list(drug2))
    }
  }
  
  # Create and return the DrugCombi object
  return(DrugCombi$new(drugs = drugs_list))
}
