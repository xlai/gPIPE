DoseConfiguration <- setRefClass("DoseConfiguration",
    fields = list(
        drugCombi = "DrugCombination",  # The DrugDose object
        currentConfig = "list",  # Current configuration of dose levels for each drug
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
            currentConfig <<- list()
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
                return(checkMonotonicity(currentConfig, increasing = isNonDecreasing))
            }
        },        
        updateConfiguration = function(newConfig, updateFunction = NULL) {
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

checkMonotonicity <- function(gamma, increasing = TRUE) {
    # Determine the number of dimensions
    numDims <- max(sapply(gamma, function(x) length(x$levels)))
    
    # Find the maximum values for each dimension
    maxValues <- sapply(1:numDims, function(d) {
        max(sapply(gamma, function(x) {
            if (length(x$levels) >= d) x$levels[d] else 0
        }))
    })
    
    # Create an array to store the gamma values, initialised with 0
    gammaArray <- array(0, dim = maxValues)
    
    # Fill the array with the corresponding validity values
    for (name in names(gamma)) {
        item <- gamma[[name]]
        indices <- item$id
        validity <- ifelse(item$validity, 1, 0) # Convert TRUE/FALSE to 1/0
        gammaArray[indices] <- validity
    }
    
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
        for (index in 1:(maxValues[d] - 1)) {
            if (compare(slices[[d]][[index]], slices[[d]][[index + 1]])) {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}