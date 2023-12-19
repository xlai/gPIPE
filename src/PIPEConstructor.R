PipeEstimator <- setRefClass("PipeEstimator",
    fields = list(
        drugCombi = "DrugCombi",  # The DrugCombi object
        trial = "PatientDataModel", # The Trial object
        validConfigs = "list",  # List of valid configurations of dose levels
        bestConfigs = "numeric" # one of the valid configurations
        epsilon = "numeric" # parameter indicating the mixing rate
        weight = "numeric" # parameter set for each dose level
    ),
    methods = list(
        initialize = function(drugCombiObject = NULL, validConfigs = NULL, epsilon = 0.5) {
            if (is.null(drugCombiObject)) {
                # Handle the case where drugDoseObject is NULL
                # This can be set to a default DrugDose object or keep it as NULL
                # Example: drugDose <<- DrugDose$new() or drugDose <<- NULL
                drugCombi <<- NULL  # or any other default initialization
                validConfigs <<- NULL
            } else {
                drugCombi <<- drugCombiObject
                validConfigs <<- validConfigs
                epsilon <<- epsilon
            }
            # Initialize the current configuration with NULL for each drug
            bestConfigs <<- rep(0, drugCombiObject$getNumberOfDoseLevels(combined = TRUE))
            weight = 1

        },    
        updatePipeEstimator <- function(p_posterior){
            log_gain_list <- lapply(
                validConfigs, 
                function(validConfiguration) {loglikelihood(p_posterior, epsilon, weight, validConfiguration)}
            )
            bestConfigs <<- validConfigs[[which.max(log_ll_list)]]
        },
        setEpsilon <- function(epsilonNew){
            epsilon <<- epsilonNew
        }
        updateWeight <- function(){
            # write a function that update the weight according to patient numbers
        }
    )

# Function to create a new PipeEstimator object
createPipeEstimator <- function(drugDoseObject = NULL, validConfigs = NULL, epsilon = 0.5) {
    return(PipeEstimator$new(drugDoseObject, validConfigs, epsilon))
}

logGain <- function(p_posterior, epsilon, weight, validConfiguration){
    log_gain_vector <- weight*validConfiguration*(log(epsilon) + log(p_posterior)) +
        weight*(1 - validConfiguration)*(log(1 - epsilon) + log(1 - p_posterior))

return(sum(logLL_vector))
}