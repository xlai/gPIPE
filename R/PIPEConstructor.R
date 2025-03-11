#' PIPE Estimator Class
#'
#' A reference class that implements the Product of Independent Beta Probabilities Escalation (PIPE)
#' design with epsilon-tapering for dose-finding studies. This class calculates posterior gain 
#' functions and determines the optimal dose configuration based on the PIPE methodology.
#'
#' @field validConfigs List of valid dose configurations
#' @field bestConfigs The optimal dose configuration (DoseConfiguration object)
#' @field validConfigs_posterior List of posterior estimates for all valid configurations
#' @field epsilon Numeric value representing the mixing rate parameter
#' @field epsilonTarget Numeric target epsilon value after tapering
#' @field weight Numeric parameter set for each dose level
#' @field taper_type Character string specifying the tapering function type
#' @field custom_taper Custom tapering function (or NULL)
#'
#' @importFrom methods setRefClass
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs ggtitle
#' @seealso \code{\link{DrugCombination-class}} for the drug combination class
#' @export
PipeEstimator <- setRefClass("PipeEstimator",
    fields = list(
        validConfigs = "list",  # List of valid configurations of dose levels
        bestConfigs = "DoseConfiguration", # one of the valid configurations
        validConfigs_posterior = "list", # posterior estimate for all valid configs
        epsilon = "numeric", # parameter indicating the mixing rate
        epsilonTarget = "numeric", #target epsilon value after tapering
        weight = "numeric", # parameter set for each dose level
        taper_type = "character",
        custom_taper = "ANY"  # Could be NULL or a function
    ),
    methods = list(
        #' @description
        #' Initialize a new PipeEstimator object
        #'
        #' @param validConfigs List of valid dose configurations
        #' @param epsilonTarget Target epsilon value after tapering
        #' @param taper_type Type of tapering function ('linear', 'quadratic', or NULL)
        #' @param custom_taper Custom tapering function (or NULL)
        initialize = function(validConfigs = NULL,
                            epsilonTarget = 0.5,
                            taper_type = NULL,
                            custom_taper = NULL) {
            if (is.null(validConfigs)) {
                validConfigs <<- NULL
                weight <<- NULL
            } else {
                validConfigs <<- validConfigs
            }
            
            # Store tapering parameters
            taper_type <<- taper_type
            custom_taper <<- custom_taper
            weight <<- 1
            
            # Set initial epsilon
            epsilonTarget <<- epsilonTarget # Store target epsilon
            epsilon <<- epsilonTarget
            if (!is.null(taper_type) || !is.null(custom_taper)) {
                # If tapering is specified, initialize with n_current = 0
                updateEpsilonWithTapering(n_total = 100, n_current = 0)
            }
        },
        
        #' @description
        #' Update the PIPE estimator with posterior probabilities
        #'
        #' @param p_posterior Vector of posterior probabilities
        #' @return List of log gain values for each valid configuration
        updatePipeEstimator = function(p_posterior){
            log_gain_list <- lapply(
                validConfigs, 
                function(validConfiguration) {logGain(p_posterior, epsilon, weight, validConfiguration$currentConfig)}
            )
            validConfigs_posterior <<- log_gain_list
            bestConfigs <<- validConfigs[[which.max(log_gain_list)]]

            return(log_gain_list)
        },
        
        #' @description
        #' Update the weight parameter
        #'
        #' @param weightNew New weight value
        updateWeight = function(weightNew){
            # write a function that update the weight according to patient numbers
            weight <<- weightNew
        },
        
        #' @description
        #' Set the epsilon parameter
        #'
        #' @param epsilonNew New epsilon value
        setEpsilon = function(epsilonNew){
            epsilon <<- epsilonNew
        },
        
        #' @description
        #' Set the target epsilon parameter
        #'
        #' @param epsilonNew New target epsilon value
        setEpsilonTarget = function(epsilonNew){
            epsilonTarget <<- epsilonNew
        },
        
        #' @description
        #' Update epsilon using the specified tapering method
        #'
        #' @param n_total Total sample size
        #' @param n_current Current sample size
        #' @return Updated epsilon value
        updateEpsilonWithTapering = function(n_total, n_current) {
            if (is.null(taper_type) && is.null(custom_taper)) {
                return(epsilon)
            }
            
            if (!is.null(custom_taper)) {
                epsilon <<- custom_taper(n_total, n_current, epsilonTarget)
            } else {
                if (taper_type == "linear") {
                    epsilon <<- (n_total - n_current) / (n_total + 1) * 0.5 + 
                               (n_current + 1) / (n_total + 1) * epsilonTarget
                } else if (taper_type == "quadratic") {
                    alpha <- n_current / n_total
                    epsilon <<- (1 - alpha)^2 * 0.5 + alpha^2 * epsilonTarget
                } else {
                    stop("Unsupported taper_type. Please use 'linear', 'quadratic', or provide a custom taper function.")
                }
            }
            return(epsilon)
        },
        
        #' @description
        #' Generate a visualization of the dose configurations
        #'
        #' @param p_posterior Vector of posterior probabilities
        #' @return A ggplot object
        plot = function(p_posterior = NULL){
            acceptable_doses <- ifelse(bestConfigs$currentConfig == 0, 'Yes', 'No')
            dose_levels <- names(bestConfigs$drugCombi$getDoseCombinationsLevel())
            # Split the vector by "." 
            # this currently only works for two drugs
            split_vector <- strsplit(dose_levels, ".", fixed = TRUE)
            # Convert the list to a matrix or dataframe
            split_df <- do.call(rbind, lapply(split_vector, function(x) {
            data.frame(drug1 = x[1], drug2 = x[2], stringsAsFactors = FALSE)
            }))

            df <- cbind(split_df, acceptable_doses, p_posterior) %>% dplyr::mutate(
                drug1 = factor(drug1, levels=drug1[1:6]),
                acceptable_doses = as.factor(acceptable_doses)
                )

            # Create the plot
            ggplot2::ggplot(df, ggplot2::aes(x = drug2, y = drug1, color = acceptable_doses, size = p_posterior)) +
            ggplot2::geom_point(alpha = 0.7, ggplot2::aes(shape = acceptable_doses))+  # Adjust alpha for point transparency if desired
            ggplot2::theme_minimal() +
            ggplot2::labs(x = "Drug 1 Dose", y = "Drug 2 Dose") +
            ggplot2::ggtitle("Joint Posterior Probability of Two Drugs")
        }
    )
)

#' Calculate log gain for PIPE method
#'
#' This function calculates the log gain for a given configuration based on
#' posterior probabilities, epsilon, and weights according to the PIPE methodology.
#'
#' @param p_posterior Vector of posterior probabilities
#' @param epsilon The epsilon parameter (mixing rate)
#' @param weight Vector of weights for each dose level
#' @param validConfiguration Binary vector representing a valid dose configuration
#'
#' @return Numeric value of the log gain
#' @keywords internal
logGain <- function(p_posterior, epsilon, weight, validConfiguration){
    log_gain_vector <- weight*validConfiguration*(log(epsilon) + log(p_posterior)) +
        weight*(1 - validConfiguration)*(log(1 - epsilon) + log(1 - p_posterior))

    return(sum(log_gain_vector))
}

#' Create a new PipeEstimator object
#'
#' This function creates a new PipeEstimator object, which implements the PIPE
#' methodology for dose-finding studies.
#'
#' @param validConfigs List of valid dose configurations
#' @param epsilonTarget Target epsilon value after tapering
#' @param taper_type Type of tapering function ('linear', 'quadratic', or NULL)
#' @param custom_taper Custom tapering function (or NULL)
#'
#' @return A PipeEstimator object
#' @export
createPipeEstimator <- function(validConfigs = NULL,
                              epsilonTarget = 0.5,
                              taper_type = NULL,
                              custom_taper = NULL) {
    return(PipeEstimator$new(
        validConfigs = validConfigs,
        epsilonTarget = epsilonTarget,
        taper_type = taper_type,
        custom_taper = custom_taper
    ))
}