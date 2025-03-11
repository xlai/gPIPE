#' Calibrate Epsilon for PIPE Estimator
#'
#' This function performs a calibration process to find the optimal epsilon value
#' for the Product of Independent Beta Probabilities Escalation (PIPE) design.
#' It evaluates different epsilon values across multiple simulated trials and selects
#' the value that maximizes the true positive rate while maintaining the false discovery
#' rate below a specified threshold.
#'
#' @param dose_configs_valid List of valid dose configurations
#' @param epsilon_range Numeric vector of epsilon values to test
#' @param num_simulations Integer number of simulations to run for each epsilon value
#' @param drugCombinationModel DrugCombinationModel object containing the model parameters
#' @param drugcombi_new DrugCombi object representing the drug combinations
#' @param patientDataModel PatientDataModel object for tracking patient data
#' @param fdr_threshold Numeric value for the maximum acceptable false discovery rate (default: 0.05)
#' @param delta Numeric value representing the toxicity increment for non-MTD doses (default: 0.2)
#' @param user_defined_scenario Optional user-defined scenario to test instead of dose_configs_valid
#' @param return_fdr_tpr Logical indicating whether to return FDR and TPR results (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item optimal_epsilon: The optimal epsilon value that satisfies the FDR threshold
#'     \item fdr_results: Vector of false discovery rates for each tested epsilon value
#'     \item tpr_results: Vector of true positive rates for each tested epsilon value
#'   }
#'
#' @details
#' The function tests each epsilon value in `epsilon_range` by running `num_simulations` 
#' trial simulations. For each simulation, it:
#'
#' 1. Creates an adjusted probability vector based on the true dose configuration
#' 2. Runs a simulated trial with the current epsilon value
#' 3. Calculates the false discovery rate (FDR) and true positive rate (TPR)
#'
#' The optimal epsilon is selected as the highest value that keeps the FDR below the 
#' specified threshold. This balances between being conservative (lower epsilon values) 
#' and efficient (higher true positive rates).
#'
#' @examples
#' \dontrun{
#' # Create required objects
#' drug_combo_model <- createDrugCombinationModel("config.yaml")
#' drug_combi <- createDrugCombi(list(Drug$new(), Drug$new()))
#' patient_data <- createPatientDataModel(drug_combi)
#' valid_configs <- list(DoseConfiguration$new(drug_combi), DoseConfiguration$new(drug_combi))
#'
#' # Calibrate epsilon
#' calibration_results <- calibrate_epsilon(
#'   dose_configs_valid = valid_configs,
#'   epsilon_range = seq(0.1, 0.9, by = 0.1),
#'   num_simulations = 100,
#'   drugCombinationModel = drug_combo_model,
#'   drugcombi_new = drug_combi,
#'   patientDataModel = patient_data
#' )
#'
#' # Access the optimal epsilon
#' optimal_epsilon <- calibration_results$optimal_epsilon
#' }
#'
#' @seealso \code{\link{PipeEstimator}} for the class that uses epsilon in the PIPE methodology
#' @import ggplot2
#' @export
calibrate_epsilon <- function(dose_configs_valid,
                              epsilon_range,
                              num_simulations,
                              drugCombinationModel,
                              drugcombi_new,
                              patientDataModel,
                              fdr_threshold = 0.05,
                              delta = 0.2,
                              user_defined_scenario = NULL,
                              return_fdr_tpr = TRUE) {
  
  # Use user-defined scenario if provided, otherwise default to dose_configs_valid
  scenarios <- if (!is.null(user_defined_scenario)) list(user_defined_scenario) else dose_configs_valid
  dose_combinations <- names(drugcombi_new$getDoseCombinationsLevel())

  # Function to calculate metrics for one simulation result
  calculate_metrics <- function(sim_result, currentConfig) {
      if (is.null(sim_result)) {
          return(list(fdr = 1, tpr = 0))
      }
      
      MTD <- sim_result$MTD
      if (all(is.na(MTD))) {
          return(list(fdr = 0, tpr = 0))
      }
      
      if (sum(currentConfig)== 0){
          return(list(fdr =  sum(currentConfig[MTD] == 0) / length(MTD), tpr = 1))
      }
      fdr <- sum(currentConfig[MTD] == 0) / length(MTD)
      tpr <- sum(currentConfig * tail(sim_result$simulation_results)[[1]]$best_config) / sum(currentConfig)
      return(list(fdr = fdr, tpr = tpr))
  }
  
  # Results storage
  results <- vector("list", length(epsilon_range))
  names(results) <- as.character(epsilon_range)

  for (epsilon in epsilon_range) {
    cat(sprintf("Testing epsilon: %.2f\n", epsilon))
    
    metrics_all_scenarios <- vector("list", length(scenarios))
    pipe_hat <- PipeEstimator$new(dose_configs_valid, epsilonTarget = epsilon, taper_type = 'linear')

    for (s in seq_along(scenarios)) {
      scenario <- scenarios[[s]]
      # Adjust prob_true based on currentConfig
      currentConfig <- scenario$currentConfig
      prob_true_adj <- currentConfig*drugCombinationModel$theta + (1 - currentConfig)*(drugCombinationModel$theta + delta)
      prob_true_list_adj <- setNames(as.list(prob_true_adj), dose_combinations)
      
      scenario_metrics <- replicate(num_simulations, {
          sim_result <- runTrialSimulation(prob_true_list_adj, drugCombinationModel, drugcombi_new, pipe_hat, patientDataModel)
          calculate_metrics(sim_result, currentConfig)
      }, simplify = FALSE)
      
      metrics_all_scenarios[[s]] <- scenario_metrics
    }      
        
    # Calculate averages
    all_fdrs <- unlist(lapply(metrics_all_scenarios, function(x) sapply(x, `[[`, "fdr")))
    all_tprs <- unlist(lapply(metrics_all_scenarios, function(x) sapply(x, `[[`, "tpr")))
    
    results[[as.character(epsilon)]] <- list(
        fdr = mean(all_fdrs, na.rm = TRUE),
        tpr = mean(all_tprs, na.rm = TRUE)
    )
    
    cat(sprintf("Epsilon: %.2f, Avg FDR: %.4f, Avg TPR: %.4f\n", 
                epsilon, results[[as.character(epsilon)]]$fdr, 
                results[[as.character(epsilon)]]$tpr))
    }
  
  # Find optimal epsilon
  valid_epsilons <- epsilon_range[sapply(results, function(x) x$fdr <= fdr_threshold)]
  optimal_epsilon <- if (length(valid_epsilons) > 0) max(valid_epsilons) else NULL
  
  if (!is.null(optimal_epsilon)) {
      cat(sprintf("Optimal epsilon found: %.2f with FDR <= %.2f\n", optimal_epsilon, fdr_threshold))
  } else {
      cat("No epsilon satisfies the FDR <= fdr_threshold criterion.\n")
  }
  
  return(list(
      optimal_epsilon = optimal_epsilon,
      fdr_results = sapply(results, `[[`, "fdr"),
      tpr_results = sapply(results, `[[`, "tpr")
  ))
}