runTrialSimulation <- function(prob_true_list, drugCombinationModel, drugcombi_new, pipe_hat, patientDataModel, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)

  simulation_results <- list()
  patientDataModel$resetTrial()
  
  while (!patientDataModel$isStoppingCriteriaMet()) {
    # Generate patient data for the current cohort at the current dose level
    outcome_prob <- as.numeric(prob_true_list[patientDataModel$currentDoseLevel])
    patientDataModel$generateRandomPatientData(outcome_prob)
    
    # Update the model based on the new patient data
    p_posterior <- drugCombinationModel$updateModel(patientDataModel)
    
    # Update PIPE estimator
    pipe_hat$updateEpsilonWithTapering(
      n_total = patientDataModel$maxCohorts,
      n_current = patientDataModel$currentCohort
    )
    temp <- pipe_hat$updatePipeEstimator(p_posterior)
    
    # Store iteration results
    simulation_results[[patientDataModel$currentCohort-1]] <- list(
        dose_level = patientDataModel$currentDoseLevel,
        p_posterior = p_posterior,
        best_config = pipe_hat$bestConfigs$currentConfig
    )
    
    current_dose_level_numeric <- drugcombi_new$getDoseCombinationsLevel(current_dose_level)
    next_dose_level_numeric <- patientDataModel$getNextDoseLevel(current_dose_level_numeric, pipe_hat, drugCombinationModel)
    if(is.na(next_dose_level_numeric)){
      break # Break if no doses found to continue the trial
    }
    next_dose_level <- names(drugcombi_new$getDoseCombinationsLevel(next_dose_level_numeric))
    patientDataModel$setCurrentDoseLevel(next_dose_level)
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