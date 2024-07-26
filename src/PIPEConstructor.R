PipeEstimator <- setRefClass("PipeEstimator",
    fields = list(
#       trial = "PatientDataModel", # The Trial object
        validConfigs = "list",  # List of valid configurations of dose levels
        bestConfigs = "DoseConfiguration", # one of the valid configurations
        validConfigs_posterior = "list", # posterior estimate for all valid configs
        epsilon = "numeric", # parameter indicating the mixing rate
        weight = "numeric" # parameter set for each dose level
    ),
    methods = list(
        initialize = function(validConfigs = NULL, epsilon = 0.5) {
            if (is.null(validConfigs)) {
                # Handle the case where drugCombiObject is NULL
                # This can be set to a default DrugDose object or keep it as NULL
                # Example: drugDose <<- DrugDose$new() or drugDose <<- NULL
#                drugCombi <<- NULL  # or any other default initialization
                validConfigs <<- NULL
                weight <<- NULL
            } else {
                validConfigs <<- validConfigs
                # Initialize the current configuration with NULL for each drug
            }
            epsilon <<- epsilon
            weight <<- 1
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
    log_gain_vector <- weight*validConfiguration*(log(epsilon) + log(1 - p_posterior)) +
        weight*(1 - validConfiguration)*(log(1 - epsilon) + log(p_posterior))

return(sum(log_gain_vector))
}