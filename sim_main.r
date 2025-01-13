# Step 1: Read in the CSV file
prob_true <- read.csv('tests/dose_levels_probability_riviere.csv')

# Step 2: Create a function to generate the prob_true_list for a given scenario
create_prob_true_list <- function(scenario_data) {
  prob_true_list <- setNames(as.list(scenario_data$probability), scenario_data$dose_level)
  return(prob_true_list)
}

# Step 3: Generate a list of prob_true_list for each scenario
prob_true_lists <- prob_true %>%
  group_by(scenario) %>%
  group_map(~ create_prob_true_list(.x))


run_simulation <- function(
  prob_true_list,
  epsilon_target,
  drugCombinationModel,
  drugcombi_new,
  pipe_hat,
  patientDataModel,
  num_simulations = 1000
) {
      # Initialize a list to store simulation results
  simulation_results <- list()
  pipe_hat$epsilonTarget <- epsilon_target

  for (i in 1:num_simulations) {
    simulation_results[[i]] <- tryCatch({
      runTrialSimulation(
        prob_true_list, drugCombinationModel, drugcombi_new,
        pipe_hat, patientDataModel
        )
    }, error = function(e) {
      if (grepl("No admissible doses found", e$message)) {
        cat(sprintf("Simulation %d stopped early: %s\n", i, e$message))
        return(NULL)
      } else {
        stop(e)
      }
    })
  }
  return(simulation_results)
}


# Step 4: Run the simulation for each scenario and store the results in a list
run_simulation_for_all_scenarios <- function(prob_true_lists, num_simulations, epsilon_target, ...) {
  results <- list()
  pipe_hat$setEpsilonTarget(epsilon_target)
  for (i in seq_along(prob_true_lists)) {
    prob_true_list <- prob_true_lists[[i]]
    results[[paste0("scenario_", i)]] <- run_simulation(
      epsilon_target = epsilon_target,          # adjust as needed
      num_simulations = num_simulations,
      prob_true_list = prob_true_list,
      drugCombinationModel,
      drugcombi_new,
      pipe_hat,
      patientDataModel
    )
  }
  return(results)
}

# Example of running the simulations with 100 simulations per scenario
results_eps0.6<- run_simulation_for_all_scenarios(prob_true_lists[1], num_simulations = 1000, epsilon_target = 0.6)

# Step 5: The results variable will now contain the simulation results for each scenario