runTrialSimulation <- function(starting_dose_level, cohort_size, max_cohorts, prob_true_list, drugCombinationModel, drugcombi_new, pipe_hat, patientDataModel, epsilon = 0.05, taper_type = 'quadratic') {
  # Initialise variables
  current_dose_level <- starting_dose_level
  simulation_results <- list()
  cohort_count <- 0  # Initialize cohort count
  
  # Create patient data model
  # patientDataModel <- createPatientDataModel(drugcombi_new, admissible_rule_list, selection_rule_list)
  
  repeat {
    # Generate patient data for the current cohort at the current dose level
    outcome_prob <- as.numeric(prob_true_list[current_dose_level])
    patientDataModel$generateRandomPatientData(current_dose_level, cohort_size, outcomeProb = outcome_prob)
    
    # Update the model based on the new patient data
    p_posterior <- drugCombinationModel$updateModel(patientDataModel)
    
    # Update PIPE estimator
    epsilon_n <- epsilon_tapering(n_total = cohort_size * max_cohorts, n_current = cohort_size * cohort_count, epsilon, taper_type = taper_type)
    pipe_hat$setEpsilon(epsilon_n)
    
    temp <- pipe_hat$updatePipeEstimator(p_posterior)
    
    # Store iteration results
    simulation_results[[length(simulation_results) + 1]] <- list(
      dose_level = current_dose_level,
      p_posterior = p_posterior,
      best_config = pipe_hat$bestConfigs$currentConfig
    )
    
    # Increment cohort count
    cohort_count <- cohort_count + 1
    
    # Check for next dose level or if maximum number of cohorts reached
    if (cohort_count >= max_cohorts) {
      break # Break if maximum number of cohorts have been recruited
    }
    
    current_dose_level_numeric <- drugcombi_new$getDoseCombinationsLevel(current_dose_level)
    next_dose_level_numeric <- patientDataModel$getNextDoseLevel(current_dose_level_numeric, pipe_hat, drugCombinationModel)
    if(is.na(next_dose_level_numeric)){
      break # Break if no doses found to continue the trial
    }
    next_dose_level <- names(drugcombi_new$getDoseCombinationsLevel(next_dose_level_numeric))
    
    current_dose_level <- next_dose_level
  }
  summStats <- patientDataModel$getSummaryStats(includeAllCombi = TRUE)
  p_posterior_mode <- 
      drugCombinationModel$calculatePosterior(drugCombinationModel$theta + 0.1, summStats) - 
      drugCombinationModel$calculatePosterior(drugCombinationModel$theta - 0.1, summStats)
  # Capture final patient data and RP2D
  final_patient_data <- patientDataModel$patientData
  tried_doses <- sort(drugcombi_new$getDoseCombinationsLevel(final_patient_data %>% pull(doseCombination) %>% unlist() %>% unique()))
  RP2D <- tried_doses[which.max(p_posterior_mode[tried_doses])]
  MTD <- patientDataModel$getMTD(pipe_hat$bestConfigs$currentConfig, drugcombi_new)
  
  return(list(
    simulation_results = simulation_results,
    final_patient_data = final_patient_data,
    RP2D = RP2D,
    MTD = MTD,
    p_posterior_mode = p_posterior_mode
  ))
}
