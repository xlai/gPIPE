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
        updateModel = function(patientDataModel) {
            # Extract summary statistics from patientDataModel
            summaryStats <- patientDataModel$getSummaryStats(includeAllCombi = TRUE)
            if (length(priorA)> 1){
                priorA_list <- setNames(as.list(priorA),names(summaryStats))
                priorB_list <- setNames(as.list(priorB),names(summaryStats))
            }
            else {
                priorA_list <- setNames(as.list(rep(priorA, length(summaryStats))),names(summaryStats))
                priorB_list <- setNames(as.list(rep(priorB, length(summaryStats))),names(summaryStats))                
            }
            if (family == "binomial") {
                updatedStats <- list()
                for (all_dose in names(summaryStats)){
                    numPatients <- summaryStats[[all_dose]]$numPatients
                    numAdverseEvents <- summaryStats[[all_dose]]$numAdverseEvents

                    posteriorA <-  priorA_list[[all_dose]] + numAdverseEvents
                    posteriorB <-  priorB_list[[all_dose]] + numPatients - numAdverseEvents
                    updatedStats[[all_dose]]$probTheta <- pbeta(theta, posteriorA, posteriorB)
                }
                relevant_names <- names(summaryStats)
                # Re-order the updatedModel_posterior list
                relevant_posterior <- updatedStats[relevant_names]
                # Assuming each element in relevant_posterior is a list that contains probTheta, extract them
                p_posterior <<- sapply(relevant_posterior, function(x) x$probTheta)

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
