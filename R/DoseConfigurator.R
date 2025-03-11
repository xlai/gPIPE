#' Dose Configuration Class
#'
#' A reference class that represents and manages dose configurations for dose-finding studies.
#' This class handles the binary representation of dose combinations, where each configuration
#' indicates whether a dose combination is acceptable (1) or unacceptable (0) based on toxicity.
#' It ensures that dose configurations satisfy monotonicity constraints.
#'
#' @field drugCombi An object of class \code{\link{DrugCombination-class}}
#' @field currentConfig A numeric vector representing the current configuration of dose levels
#'        (1 for acceptable doses, 0 for unacceptable doses)
#' @field isNonDecreasing Logical indicating whether the configuration follows a non-decreasing pattern
#'        (TRUE) or non-increasing pattern (FALSE)
#' @field isValid Logical indicating whether the current configuration satisfies monotonicity constraints
#'
#' @seealso \code{\link{PipeEstimator-class}} for the class that uses dose configurations to implement
#'          the PIPE methodology, \code{\link{DrugCombination-class}} for the underlying drug combination class
#' @importFrom methods setRefClass
#' @export
DoseConfiguration <- setRefClass("DoseConfiguration",
    fields = list(
        drugCombi = "DrugCombination",  # The DrugDose object
        currentConfig = "numeric",  # Current configuration of dose levels for each drug
        isNonDecreasing = "logical",  # TRUE for non-decreasing, FALSE for non-increasing
        isValid = "logical" # TRUE for valid dose config, FALSE otherwise
    ),
    methods = list(
        #' @description
        #' Initialize a new DoseConfiguration object
        #'
        #' @param drugCombiObject An object of class DrugCombination
        #' @param isNonDecreasing Logical indicating monotonicity direction (TRUE for non-decreasing)
        #' @param currentConfig Optional numeric vector representing the initial configuration
        initialize = function(drugCombiObject = NULL, isNonDecreasing = TRUE, currentConfig = NULL) {
            if (is.null(drugCombiObject)) {
                # Handle the case where drugDoseObject is NULL
                drugCombi <<- DrugCombi$new()  # or any other default initialization
                currentConfig <<- numeric(length = length(drugCombi$doseCombinations))
            } else {
                drugCombi <<- drugCombiObject
                currentConfig <<- currentConfig
            }
            isNonDecreasing <<- isNonDecreasing
            isValid <<- .self$isValidConfiguration()
        },
        
        #' @description
        #' Check if the current configuration satisfies monotonicity constraints
        #'
        #' @return Logical indicating whether the configuration is valid
        isValidConfiguration = function() {
            # Check the validity of its current configurations and return logical value
            if (length(currentConfig) == 0) {
                return(FALSE)  # Invalid configuration length
            }
            else{
                return(checkMonotonicity(currentConfig, drugCombi, increasing = isNonDecreasing))
            }
        },
        
        #' @description
        #' Update the current configuration
        #'
        #' @param newConfig Numeric vector representing the new configuration
        updateConfiguration = function(newConfig) {
                # Manually set the new configuration
                currentConfig <<- newConfig
                isValid <<- .self$isValidConfiguration()
        },
        
        #' @description
        #' Get the number of dose levels in the configuration
        #'
        #' @return Numeric value representing the number of dose levels
        getNumberOfDoseLevels = function() {
            return(length(currentConfig))
        },
        
        #' @description
        #' Get indices of acceptable doses (those with value 1)
        #'
        #' @return Numeric vector of indices for acceptable doses
        getAcceptableDoses = function() {
            return(which(currentConfig == 1))
        },
        
        #' @description
        #' Get indices of unacceptable doses (those with value 0)
        #'
        #' @return Numeric vector of indices for unacceptable doses
        getUnacceptableDoses = function() {
            return(which(currentConfig == 0))
        },
        
        #' @description
        #' Generate a string representation of the configuration
        #'
        #' @return Character string representing the configuration
        toString = function() {
            return(paste(currentConfig, collapse = ""))
        },
        
        #' @description
        #' Print the dose configuration
        print = function() {
            cat("Dose Configuration:\n")
            cat("  Valid:", ifelse(isValid, "Yes", "No"), "\n")
            cat("  Monotonicity:", ifelse(isNonDecreasing, "Non-decreasing", "Non-increasing"), "\n")
            
            # Get drug combination names
            doseNames <- names(drugCombi$getDoseCombinationsLevel())
            
            # Create a data frame for pretty printing
            config_df <- data.frame(
                Dose = doseNames,
                Acceptable = ifelse(currentConfig == 1, "Yes", "No")
            )
            
            # Print the data frame
            print(config_df)
        }
    )
)

#' Slice an array by each dimension
#'
#' This helper function slices a multi-dimensional array along each dimension.
#'
#' @param array The array to slice
#'
#' @return A list of lists, where each inner list contains slices along a dimension
#' @keywords internal
slice_array_by_dimension <- function(array) {
  dimensions <- dim(array)
  num_dimensions <- length(dimensions)

  # Use lapply to iterate over dimensions
  all_slices <- lapply(1:num_dimensions, function(dim) {
    # Use lapply to generate slices for each index in the dimension
    lapply(1:dimensions[dim], function(i) {
      index_expr <- lapply(1:num_dimensions, function(x) ifelse(x == dim, i, TRUE))
      do.call("[", c(list(array), index_expr))
    })
  })

  return(all_slices)
}

#' Check monotonicity of a dose configuration
#'
#' This function verifies that a dose configuration satisfies monotonicity constraints.
#' For non-decreasing configurations, higher doses should not be acceptable if lower doses are unacceptable.
#' For non-increasing configurations, lower doses should not be acceptable if higher doses are unacceptable.
#'
#' @param currentConfig Numeric vector representing the configuration to check
#' @param drugCombi DrugCombination object
#' @param increasing Logical indicating whether to check for non-decreasing (TRUE) or non-increasing (FALSE) monotonicity
#'
#' @return Logical indicating whether the configuration satisfies monotonicity constraints
#' @keywords internal
checkMonotonicity <- function(currentConfig, drugCombi, increasing = TRUE) {
    # Determine the number of dimensions
    numDoseLevels <- drugCombi$getNumberOfDoseLevels(combined = FALSE)
    numDims <- length(numDoseLevels)
    indices <- drugCombi$getDoseCombinationsLevel()
    # Create an array to store the gamma values, initialised with 0
    gammaArray <- matrix(currentConfig, ncol = numDoseLevels[1], byrow = FALSE)
    
    # Pre-compute slices for each dimension
    slices <- slice_array_by_dimension(gammaArray)
    
    # Define the comparison function based on the type of monotonicity
    compare <- if (increasing) {
        function(slice1, slice2) { any(slice2 < slice1, na.rm = TRUE) }
    } else {
        function(slice1, slice2) { any(slice2 > slice1, na.rm = TRUE) }
    }

    # Check monotonicity in each dimension
    for (d in 1:numDims) {
        for (index in 1:(numDoseLevels[-d] - 1)) {
            if (compare(slices[[d]][[index]], slices[[d]][[index + 1]])) {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}

#' Create a new DoseConfiguration object
#'
#' This function creates a new DoseConfiguration object, which represents and manages
#' dose configurations in dose-finding studies.
#'
#' @param drugCombiObject An object of class DrugCombination
#' @param isNonDecreasing Logical indicating monotonicity direction (TRUE for non-decreasing)
#' @param currentConfig Optional numeric vector representing the initial configuration
#'
#' @return A DoseConfiguration object
#' @seealso \code{\link{DoseConfiguration-class}} for the class documentation
#' @export
createDoseConfig <- function(drugCombiObject = NULL, isNonDecreasing = TRUE, currentConfig = NULL) {
    return(DoseConfiguration$new(drugCombiObject, isNonDecreasing, currentConfig))
}