PipeEstimator <- setRefClass("PipeEstimator",
    fields = list(
        drugDose = "DrugDose",  # The DrugDose object
        currentConfig = "list",  # Current configuration of dose levels for each drug
        isNonDecreasing = "logical"  # TRUE for non-decreasing, FALSE for non-increasing
    ),
    methods = list(
        initialize = function(drugDoseObject = NULL, isNonDecreasing = TRUE) {
            if (is.null(drugDoseObject)) {
                # Handle the case where drugDoseObject is NULL
                # This can be set to a default DrugDose object or keep it as NULL
                # Example: drugDose <<- DrugDose$new() or drugDose <<- NULL
                drugDose <<- NULL  # or any other default initialization
            } else {
                drugDose <<- drugDoseObject
            }
            isNonDecreasing <<- isNonDecreasing
            # Initialize the current configuration with NULL for each drug
            if (!is.null(drugDose)) {
                currentConfig <<- lapply(drugDose$drugs, function(x) NULL)
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
                for (drug_index in seq_along(gamma)) {
                    currentNumericOrder <- if (is.null(currentConfig[[drug_index]])) {
                        if (isNonDecreasing) 1 else max(drugDose$getNumericOrder(names(drugDose$drugs)[drug_index]))
                    } else {
                        drugDose$getNumericOrder(names(drugDose$drugs)[drug_index], currentConfig[[drug_index]])
                    }
                    
                    newNumericOrder <- drugDose$getNumericOrder(names(drugDose$drugs)[drug_index], gamma[[drug_index]])

                    if ((isNonDecreasing && newNumericOrder < currentNumericOrder) ||
                        (!isNonDecreasing && newNumericOrder > currentNumericOrder)) {
                        return(FALSE)
                    }
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

# Function to create a new PipeEstimator object
createPipeEstimator <- function(drugDoseObject = NULL, isNonDecreasing = TRUE) {
    return(PipeEstimator$new(drugDoseObject, isNonDecreasing))
}