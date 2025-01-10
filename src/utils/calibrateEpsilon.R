source('src/utils/map_vector.R')

calibrate_epsilon <- function(dose_configs_valid, 
                              epsilon_range, 
                              num_simulations, 
                              starting_dose_level, 
                              cohort_size, 
                              max_cohorts, 
                              drugCombinationModel, 
                              drugcombi_new, 
                              pipe_hat, 
                              fdr_threshold = 0.05, 
                              user_defined_scenario = NULL,
                              return_fdr_tpr = TRUE) {
  
  # Use user-defined scenario if provided, otherwise default to dose_configs_valid
  scenarios <- if (!is.null(user_defined_scenario)) list(user_defined_scenario) else dose_configs_valid
  
  # Initialize a list to store the FDR and TPR results for each epsilon
  fdr_results <- list()
  tpr_results <- list()
  
  for (epsilon in epsilon_range) {
    cat(sprintf("Testing epsilon: %.2f\n", epsilon))
    
    # Initialize variables to store FDR and TPR for all scenarios
    fdr_scenario <- numeric()
    tpr_scenario <- numeric()
    
    for (scenario in scenarios) {
      patientDataModel <- createPatientDataModel(drugcombi_new, admissible_rule_list, selection_rule_list)
      # Adjust prob_true based on currentConfig
      currentConfig <- scenario$currentConfig
      mtd_idx <- patientDataModel$getMTD(currentConfig, drugcombi_new)
      currentConfig1 <- currentConfig
#      currentConfig1[mtd_idx] <- drugCombinationModel$theta
      prob_true_adj <- map_vector2(1 - currentConfig1, drugCombinationModel$theta, drugCombinationModel$theta)
      prob_true_list_adj <- setNames(
        as.list(prob_true_adj),
        names(drugcombi_new$getDoseCombinationsLevel())
      )
      
      # Initialize a list to store simulation results
      simulation_results <- list()
      
      for (i in 1:num_simulations) {
        patientDataModel <- createPatientDataModel(drugcombi_new, admissible_rule_list, selection_rule_list)
        sim_result <- tryCatch({
          runTrialSimulation(starting_dose_level, cohort_size, max_cohorts, 
                            prob_true_list_adj, drugCombinationModel, drugcombi_new, 
                            pipe_hat, patientDataModel, epsilon = epsilon)
        }, error = function(e) {
          if (grepl("No admissible doses found", e$message)) {
            cat(sprintf("Simulation %d stopped early: %s\n", i, e$message))
            return(NULL)
          } else {
            stop(e)
          }
        })
        
        # Calculate FDR and TPR if sim_result is not NULL
        if (!is.null(sim_result)) {
          MTD <- sim_result$MTD
          
          if (all(is.na(MTD))) {  # No MTD found
            fdr <- NA
            tpr <- 0
          } else {
            fdr <- sum(currentConfig[MTD]==0) / length(MTD)
            if (sum(currentConfig) == length(currentConfig)) {
              tpr <- sum(currentConfig * tail(sim_result$simulation_results)[[1]]$best_config)/ sum(currentConfig)              
            } else{
              tpr <- sum(currentConfig * tail(sim_result$simulation_results)[[1]]$best_config)/ sum(currentConfig)
            }
          }
          
          fdr_scenario <- c(fdr_scenario, fdr)
          tpr_scenario <- c(tpr_scenario, tpr)
        } else {
          # No admissible dose found in the simulation, setting FDR = 1 and TPR = 0
          fdr_scenario <- c(fdr_scenario, 1)
          tpr_scenario <- c(tpr_scenario, 0)
        }
      }
    }
    
    # Calculate average FDR and TPR for this epsilon
    avg_fdr <- mean(fdr_scenario, na.rm = TRUE)
    avg_tpr <- mean(tpr_scenario, na.rm = TRUE)
    
    cat(sprintf("Epsilon: %.2f, Avg FDR: %.4f, Avg TPR: %.4f\n", epsilon, avg_fdr, avg_tpr))
    
    fdr_results[[as.character(epsilon)]] <- avg_fdr
    tpr_results[[as.character(epsilon)]] <- avg_tpr
  }
  
  # Identify the largest epsilon with FDR <= fdr_threshold
  valid_epsilons <- names(fdr_results)[sapply(fdr_results, function(x) x <= fdr_threshold)]
  optimal_epsilon <- if (length(valid_epsilons) > 0) {
    max(as.numeric(valid_epsilons))
  } else {
    NULL
  }
  
  if (!is.null(optimal_epsilon)) {
    cat(sprintf("Optimal epsilon found: %.2f with FDR <= %.2f\n", optimal_epsilon, fdr_threshold))
  } else {
    cat("No epsilon satisfies the FDR <= fdr_threshold criterion.\n")
  }
  
  # Return a list with optimal epsilon, FDR results, and TPR results
  return(list(optimal_epsilon = optimal_epsilon, 
              fdr_results = fdr_results, 
              tpr_results = tpr_results))
}