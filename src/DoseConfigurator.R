DoseConfiguration <- setRefClass("DoseConfiguration",
    fields = list(
        drugCombi = "DrugCombination",  # The DrugDose object
        currentConfig = "numeric",  # Current configuration of dose levels for each drug
        isNonDecreasing = "logical",  # TRUE for non-decreasing, FALSE for non-increasing
        isValid = "logical" # TRUE for valid dose config, FALSE otherwise
    ),
    methods = list(
        initialize = function(drugCombiObject = NULL, isNonDecreasing = TRUE, currentConfig = NULL) {
            if (is.null(drugCombiObject)) {
                # Handle the case where drugDoseObject is NULL
                # This can be set to a default DrugDose object or keep it as NULL
                # Example: drugDose <<- DrugDose$new() or drugDose <<- NULL
                drugCombi <<- DrugCombi$new()  # or any other default initialization
            } else {
                drugCombi <<- drugCombiObject
            }
            currentConfig <<- numeric(length = length(drugCombi$doseCombinations))
            isNonDecreasing <<- isNonDecreasing
            isValid <<- .self$isValidConfiguration()
        },
        isValidConfiguration = function() {
            # Check the validity of its current configurations and return logcial value
            # Function to check a single configuration
            if (length(currentConfig) == 0) {
                return(FALSE)  # Invalid configuration length
            }
            else{
                return(checkMonotonicity(currentConfig, drugCombi, increasing = isNonDecreasing))
            }
        },        
        updateConfiguration = function(newConfig) {
                # Manually set the new configuration
                currentConfig <<- newConfig
                isValid <<- .self$isValidConfiguration()
        }        
        # ... (other methods if necessary)
    )
)

# Function to slice an array by each dimension using lapply
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

checkMonotonicity <- function(currentConfig, drugCombi, increasing = TRUE) {
    # Determine the number of dimensions
    numDoseLevels <- drugCombi$getNumberOfDoseLevels(combined = FALSE)
    numDims <- length(numDoseLevels)
    indices <- drugCombi$getDoseCombinationsLevel()
    # Create an array to store the gamma values, initialised with 0
    gammaArray <- array(currentConfig, dim = numDoseLevels)
    
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
        for (index in 1:(numDoseLevels[d] - 1)) {
            if (compare(slices[[d]][[index]], slices[[d]][[index + 1]])) {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}