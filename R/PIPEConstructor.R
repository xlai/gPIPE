PipeEstimator <- setRefClass("PipeEstimator",
    fields = list(
#       trial = "PatientDataModel", # The Trial object
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
        updatePipeEstimator = function(p_posterior){
            log_gain_list <- lapply(
                validConfigs, 
                function(validConfiguration) {logGain(p_posterior, epsilon, weight, validConfiguration$currentConfig)}
            )
            validConfigs_posterior <<- log_gain_list
            bestConfigs <<- validConfigs[[which.max(log_gain_list)]]

            return(log_gain_list)
        },
        updateWeight = function(weightNew){
            # write a function that update the weight according to patient numbers
            weight <<- weightNew
        },        
        setEpsilon = function(epsilonNew){
            epsilon <<- epsilonNew
        },
        setEpsilonTarget = function(epsilonNew){
            epsilonTarget <<- epsilonNew
        },        
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

            df <- cbind(split_df, acceptable_doses, p_posterior) %>% mutate(
                drug1 = factor(drug1, levels=drug1[1:6]),
                acceptable_doses = as.factor(acceptable_doses)
                )

            # Create the plot
            ggplot(df, aes(x = drug2, y = drug1, color = acceptable_doses, size = p_posterior)) +
            geom_point(alpha = 0.7, aes(shape = acceptable_doses))+  # Adjust alpha for point transparency if desired
#            scale_color_gradient(low = "blue", high = "red") +  # Customize color gradient
            theme_minimal() +
            labs(x = "Drug 1 Dose", y = "Drug 2 Dose") +
            ggtitle("Joint Posterior Probability of Two Drugs")
        }

    )
)



logGain <- function(p_posterior, epsilon, weight, validConfiguration){
    log_gain_vector <- weight*validConfiguration*(log(epsilon) + log(p_posterior)) +
        weight*(1 - validConfiguration)*(log(1 - epsilon) + log(1 - p_posterior))

return(sum(log_gain_vector))
}