PipeEstimator <- setRefClass("PipeEstimator",
    fields = list(
        drugCombi = "DrugCombi",  # The DrugDose object
        currentConfig = "list",  # Current configuration of dose levels for each drug
        isNonDecreasing = "logical"  # TRUE for non-decreasing, FALSE for non-increasing
    ),
    methods = list(
        initialize = function(drugCombiObject = NULL, isNonDecreasing = TRUE) {
            if (is.null(drugCombiObject)) {
                # Handle the case where drugDoseObject is NULL
                # This can be set to a default DrugDose object or keep it as NULL
                # Example: drugDose <<- DrugDose$new() or drugDose <<- NULL
                drugCombi <<- NULL  # or any other default initialization
            } else {
                drugCombi <<- drugCombiObject
            }
            isNonDecreasing <<- isNonDecreasing
            # Initialize the current configuration with NULL for each drug
            if (!is.null(drugCombiObject)) {
                currentConfig <<- lapply(drugCombiObject$doseCombinations, function(x) NULL)
            } else {
                currentConfig <<- NULL
            }
        },
        isValidConfiguration = function(gammaConfigs) {
            # Check the validity of multiple gamma configurations
            # gammaConfigs: A list of configurations, each configuration being a list of dose levels for each drug

            # Function to check a single configuration
            checkSingleConfig <- function(gamma) {
                if (length(gamma) != length(currentConfig)) {
                    return(FALSE)  # Invalid configuration length
                }
                else{
                    return(checkMonotonicity(gamma, increasing = isNonDecreasing))
                }
                return(TRUE)
            }

            # Apply the check to each configuration in the list
            return(sapply(gammaConfigs, checkSingleConfig))
        },        
        updateConfiguration = function(newConfig, updateFunction = NULL) {
            if (!is.null(updateFunction) && is.function(updateFunction)) {
                # Update configuration using the provided function
                currentConfig <<- updateFunction(currentConfig, newConfig)
            } else {
                # Manually set the new configuration
                currentConfig <<- newConfig
            }
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

# Function to create a new PipeEstimator object
createPipeEstimator <- function(drugDoseObject = NULL, isNonDecreasing = TRUE) {
    return(PipeEstimator$new(drugDoseObject, isNonDecreasing))
}

loglikelihood <- function(p_posterior, epsilon, weight, pipe_configuration){
    logLL_vector <- weight*pipe_configuration*(log(epsilon) + log(p_posterior)) +
        weight*(1 - pipe_configuration)*(log(1 - epsilon) + log(1 - p_posterior))

return(sum(logLL_vector))
}