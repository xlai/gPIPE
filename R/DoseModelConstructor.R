#' Drug Combination Statistical Model
#'
#' A reference class that implements the statistical model behind dose combinations
#' for dose-finding studies. This class provides methods for Bayesian inference using 
#' beta priors and binomial likelihood to calculate posterior probabilities of toxicity 
#' at different dose combinations.
#'
#' @field theta Numeric value representing the target toxicity probability threshold
#' @field priorA Numeric vector of beta prior alpha parameters for each dose combination
#' @field priorB Numeric vector of beta prior beta parameters for each dose combination
#' @field family Character string specifying the distribution family ("binomial" or "gaussian")
#' @field p_posterior Numeric vector of posterior probabilities
#'
#' @importFrom methods setRefClass
#' @importFrom yaml yaml.load_file
#' @seealso \code{\link{DrugCombination-class}} for the drug combination class that this model works with
#' @export
DrugCombinationModel <- setRefClass("DrugCombinationModel",
    fields = list(
        theta = "numeric",       # Target toxicity probability threshold
        priorA = "numeric",      # Beta prior alpha parameters
        priorB = "numeric",      # Beta prior beta parameters
        family = "character",    # Distribution family
        p_posterior = "numeric"  # Posterior probabilities
    ),
    methods = list(
        #' @description
        #' Initialize a new DrugCombinationModel object
        #'
        #' @param yamlConfig Path to YAML configuration file containing model parameters
        #' @param family Character string specifying the distribution family ("binomial" or "gaussian")
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
        
        #' @description
        #' Calculate posterior probability for one or all doses
        #'
        #' @param probValue Numeric value at which to evaluate the posterior CDF
        #' @param summaryStats List with summary statistics for each dose combination
        #' @param doseName Optional character string for a specific dose name
        #' @return Numeric vector of posterior probabilities
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
        
        #' @description
        #' Update the model with patient data
        #'
        #' @param patientDataModel A obejct of class PatientDataModel containing the trial data
        #' @return Numeric vector of updated posterior probabilities
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

#' Create a new DrugCombinationModel object
#'
#' This function creates a new DrugCombinationModel object, which implements
#' the statistical model for a trial.
#'
#' @param yamlConfig Path to YAML configuration file containing model parameters
#' @param family Character string specifying the distribution family ("binomial" or "gaussian")
#'
#' @return A DrugCombinationModel object
#' @export
createDrugCombinationModel <- function(yamlConfig, family = "binomial") {
    return(DrugCombinationModel$new(yamlConfig, family))
}
