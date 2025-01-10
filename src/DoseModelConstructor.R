DrugCombinationModel <- setRefClass("DrugCombinationModel",
    fields = list(
        theta = "numeric",
        priorA = "numeric",
        priorB = "numeric",
        family = "character",
        p_posterior = "numeric"
    ),
    methods = list(
        initialize = function(yamlConfig, family = "binomial") {
            config <- yaml::yaml.load_file(yamlConfig)
            theta <<- config$theta
            
            if (!is.null(config$priorA) && !is.null(config$priorB)) {
                priorA <<- config$priorA
                priorB <<- config$priorB
            } else if (!is.null(config$priorMedian) && !is.null(config$priorStrength)) {
                # Calculate priorA and priorB using prior median and sample size
                if (config$priorStrength == 'Strong') {
                    priorSampleSize = rep(1, length(config$priorMedian))
                }
                else if (config$priorStrength == 'Weak'){
                    priorSampleSize = rep(1/length(config$priorMedian), length(config$priorMedian))
                }
                else {
                     stop("Please specify strengths of prior.")
                }
                betas <- find_beta_parameters(config$priorMedian, priorSampleSize)
                priorA <<- betas$a
                priorB <<- betas$b
            } else {
                stop("Either specify priorA and priorB directly or provide priorMedian with priorSampleSize.")
            }
            
            family <<- family
        },
        calculatePosterior = function(probValue, summaryStats, doseName = NULL) {
            # This function calculates posterior probabilities.
            # If doseName is provided, calculate for that dose. Otherwise, calculate for all doses.

            calculate_for_dose <- function(dose) {
                numPatients <- summaryStats[[dose]]$numPatients
                numAdverseEvents <- summaryStats[[dose]]$numAdverseEvents
                
                # Get the priorA and priorB values for the dose
                priorA_val <- ifelse(length(priorA) > 1, priorA[dose], priorA)
                priorB_val <- ifelse(length(priorB) > 1, priorB[dose], priorB)
                
                # Calculate posterior parameters
                posteriorA <- priorA_val + numAdverseEvents
                posteriorB <- priorB_val + numPatients - numAdverseEvents
                
                # Return the posterior probability at the specified probValue
                return(pbeta(probValue, posteriorA, posteriorB))
            }

            if (is.null(doseName)) {
                # If no doseName is provided, calculate posterior for all doses
                return(sapply(names(summaryStats), calculate_for_dose))
            } else {
                # If a specific doseName is provided, calculate for that dose
                return(calculate_for_dose(doseName))
            }
        },
        
        updateModel = function(patientDataModel) {
            # Extract summary statistics from patientDataModel
            summaryStats <- patientDataModel$getSummaryStats(includeAllCombi = TRUE)
            
            if (length(priorA) > 1){
                priorA_list <- setNames(as.list(priorA), names(summaryStats))
                priorB_list <- setNames(as.list(priorB), names(summaryStats))
            } else {
                priorA_list <- setNames(as.list(rep(priorA, length(summaryStats))), names(summaryStats))
                priorB_list <- setNames(as.list(rep(priorB, length(summaryStats))), names(summaryStats))                
            }

            if (family == "binomial") {
                # Use the calculatePosterior function for all doses
                p_posterior <<- calculatePosterior(theta, summaryStats)                
                return(p_posterior)
            } else if (family == "gaussian") {
                # Gaussian family update logic (placeholder)
                # ...
            } else {
                stop("Invalid family specified.")
            }
        }
    )
)

# Function to create a new DrugCombinationModel object
createDrugCombinationModel <- function(yamlConfig, family = "binomial") {
    return(DrugCombinationModel$new(yamlConfig, family))
}
