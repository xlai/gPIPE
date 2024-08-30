library(doFuture)

# Parameters for simulation
num_simulations <- 100
starting_dose_level <- "dose1.dose1"
cohort_size <- 2
max_cohorts <- 30

# Ensure the necessary objects and lists are defined or loaded here
patientDataModel <- createPatientDataModel(drugcombi_new, admissible_rule_list, selection_rule_list)

plan(multisession)

options(future.globals.onReference = "warning")

# Run simulations in parallel
simulation_outputs <- foreach(
  i = 1:num_simulations, 
  .combine = 'list', .options.future = list(seed = TRUE, noexport=c('patientDataModel','drugCombinationModel', 'drugcombi_new'))
  ) %dofuture% {
    source('src/DrugCombiConstructor.R')
    source('src/TrialConstructor.R')
    source('src/DoseModelConstructor.R')
    drugcombi_new <- DrugCombi$new(drugs=drug_list)     
    drugCombinationModel <- createDrugCombinationModel("./src/drugPrior_6level.yaml")
    patientDataModel <- createPatientDataModel(drugcombi_new, admissible_rule_list, selection_rule_list)
    pipe_hat <- PipeEstimator(dose_configs_valid, epsilon = 0.1)

    tryCatch({
      # Call your function with the desired parameters
      runTrialSimulation(starting_dose_level, cohort_size, max_cohorts, 
                        prob_true_list, drugCombinationModel, drugcombi_new, 
                        pipe_hat, patientDataModel)
    }, error = function(e) {
      # Handle the case where no admissible dose is found
      if (grepl("No admissible doses found", e$message)) {
        cat(sprintf("Simulation %d stopped early: %s\n", i, e$message))
        return(NULL)  # Return NULL or any other indicator of early stop
      } else {
        stop(e)  # Re-throw other errors that are not handled
      }
    })
#  runTrialSimulation(starting_dose_level, cohort_size, max_cohorts, prob_true_list, drugCombinationModel, drugcombi_new, pipe_hat, admissible_rule_list, selection_rule_list)
}

run_simulation <- function(
  dose_configs_valid, 
  epsilon_final,
  num_simulations, 
  starting_dose_level, 
  cohort_size, 
  max_cohorts, 
  drugCombinationModel, 
  drugcombi_new, 
  pipe_hat, 
  taper_type = 'quadratic'
  ) {
      # Initialize a list to store simulation results
  simulation_results <- list()

  for (i in 1:num_simulations) {
    patientDataModel <- createPatientDataModel(drugcombi_new, admissible_rule_list, selection_rule_list)
    simulation_results[[i]] <- tryCatch({
      runTrialSimulation(starting_dose_level, cohort_size, max_cohorts, 
                        prob_true_list, drugCombinationModel, drugcombi_new, 
                        pipe_hat, patientDataModel, epsilon = epsilon_final, taper_type = taper_type)
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
        
# Assuming your sim_output list is already loaded
final_patient_data <- sim_output$final_patient_data

calculate_cumulative_dlt <- function(data, drugA, drugB, cohort_limit, expand_grid = FALSE) {
  # Filter data up to the specified cohort
  filtered_data <- data %>% filter(cohort <= cohort_limit)
  
  # Calculate cumulative DLT proportion for each dose combination level up to the specified cohort
  cumulative_dlt_data <- filtered_data %>%
    arrange(cohort) %>%
    group_by(doseCombination) %>%
    summarise(
      cumulative_count = n(),
      cumulative_DLT = sum(outcome),
      cumulative_DLT_proportion = cumulative_DLT / cumulative_count
      ) %>%
    ungroup() 
  
  # Convert doseCombination to separate dose levels
  cumulative_dlt_data <- cumulative_dlt_data %>%
    tidyr::separate(doseCombination, into = c("Drug1", "Drug2"), sep = "\\.") %>%
    mutate(Drug1 = factor(Drug1, levels = names(drugA$doses)),
           Drug2 = factor(Drug2, levels = names(drugB$doses)))
  
  if (expand_grid == TRUE){
    all_combinations <- expand.grid(Drug1 = names(drugA$doses), Drug2 = names(drugB$doses))
    cumulative_dlt_data <- all_combinations %>% 
      left_join(cumulative_dlt_data, by = c("Drug1", "Drug2"))
#      mutate(cumulative_DLT = ifelse(is.na(cumulative_DLT), 0, cumulative_DLT),
#      cumulative_count = ifelse(is.na(cumulative_count), 0, cumulative_count),
#      cumulative_DLT_proportion = ifelse(is.na(cumulative_DLT_proportion), 0, cumulative_DLT_proportion))

  }
  
  return(cumulative_dlt_data)
}

# Define the function
prepare_step_data <- function(sim_output, drugA, drugB, cohort_limit) {
  # Add toxicity status to the dataframe
  pipe_estimate <- sim_output$simulation_results[[cohort_limit]]$best_config
  nlevel_drugA <- length(drugA$doses)
  nlevel_drugB <- length(drugB$doses)

  xlevel <- 1:(nlevel_drugA + 1) - 0.5
  ylevel <- c(apply(matrix(pipe_estimate, ncol = nlevel_drugA),2,function(i){min(which(i==1),nlevel_drugB + 1)-0.5}), 0.5)

  # Prepare the step data
  step_data <- data.frame(
    x = xlevel,
    y = ylevel
  )

  return(step_data)
}


# Create the ggplot

plot_dlt_proportion <- function(sim_output, drugA, drugB, cohort_seq) {
  # Calculate cumulative DLT and prepare step data
  cumulative_dlt_df <- calculate_cumulative_dlt(sim_output$final_patient_data, drugA, drugB, cohort_seq, expand_grid = TRUE)
  step_data <- prepare_step_data(sim_output, drugA, drugB, cohort_seq)
  n_doses_drugA <- as.numeric(length(drugA$doses))
  n_doses_drugB <- as.numeric(length(drugB$doses))  
  # Create the plot
  p <- ggplot() +
    geom_tile(data = cumulative_dlt_df, aes(x = Drug1, y = Drug2, fill = cumulative_DLT_proportion), alpha = 0.5, size = 2) +
    scale_x_discrete(breaks = names(drugA$doses), expand = c(0, 0), drop = FALSE) +
    scale_y_discrete(breaks = names(drugB$doses), expand = c(0, 0), drop = FALSE) +
    geom_step(aes(x = x, y = y), data = step_data, linewidth = 1.5, color = "black", linetype = 2) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
                         limits = c(0, 1),
                         na.value = 'white',
                         name = "DLT Proportion") +
#    scale_x_continuous(breaks = 1:n_doses_drugA, labels = names(drugA$doses), limits = c(0, n_doses_drugA + 1), expand = c(0, 0)) +
#    scale_y_continuous(breaks = 1:n_doses_drugB, labels = names(drugB$doses), limits = c(0, n_doses_drugB + 1), expand = c(0, 0)) +
    geom_text(data = cumulative_dlt_df %>% filter(!is.na(cumulative_DLT_proportion)), aes(x = as.numeric(Drug1), y = as.numeric(Drug2), label = paste(cumulative_count, "(", cumulative_DLT, ")", sep = "")), color = "black") +
    labs(x = "Drug 1", y = "Drug 2") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


# For animation, create a list of plots for each cohort sequence
cohort_seqs <- 1:max_cohorts # Example sequence
plots <- lapply(cohort_seqs, function(cohort_seq) {
  plot_dlt_proportion(z1[[12]], drugA, drugB, cohort_seq)
})
# Extract the legend from one of the plots
legend <- cowplot::get_legend(plots[[1]] + theme(legend.position = "bottom"))

# Remove legends from all plots
plots <- lapply(plots, function(plot) plot + theme(legend.position = "none" ))

# Combine plots into a single page using cowplot
combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 4, align = 'v')

# Add a common title
title <- cowplot::ggdraw() + cowplot::draw_label("Cumulative DLT Proportion by Dose Combination", fontface = 'bold')

# Combine title, legend, and plots
final_plot <- cowplot::plot_grid(title, combined_plot, legend, ncol = 1, rel_heights = c(0.1, 1, 0.1))

ggsave(final_plot, file = '~/Downloads/pipe_plot_eps_0.3.pdf', height = 40, width = 30)


# Combine the plots into an animation
anim <- plots[[1]] + 
  transition_states(cohort_seqs, transition_length = 2, state_length = 1) +
  enter_fade() + exit_fade()

# Save the animation
anim_save("dlt_proportion_animation.gif", anim)



# Function to create intervals around the target DLT probability
create_intervals <- function(target_DLT_prob, interval_width) {
  lower_bound <- max(0, target_DLT_prob - interval_width / 2)
  upper_bound <- min(1, target_DLT_prob + interval_width / 2)
  
  # Define the intervals
  intervals <- c(0, lower_bound - interval_width, lower_bound, upper_bound, upper_bound + interval_width, 1)
  intervals <- unique(pmin(pmax(intervals, 0), 1))
  
  return(intervals)
}

# Function to create summary table with mode switching
create_summary_table <- function(prob_true_list, target_DLT_prob, indifference_interval, patientData, mode = "MTD") {
  
  # Use the provided create_intervals function to generate intervals
  intervals <- create_intervals(target_DLT_prob, indifference_interval)
  interval_labels <- sapply(seq_along(intervals[-1]), function(i) {
    paste0("(", intervals[i], ",", intervals[i+1], "]")
  })
  
  # Define the lower and upper bounds for the MTD based on the first interval around the target DLT probability
  lower_bound <- intervals[3]
  upper_bound <- intervals[4]
  
  # Map dose combinations to their true DLT probabilities to a data frame
  dose_probs_df <- data.frame(
    doseCombination = names(prob_true_list),
    DLT_prob = as.numeric(prob_true_list),
    stringsAsFactors = FALSE
  )
  
  # Add true DLT probabilities to patientData
  patientData <- patientData %>%
    mutate(doseCombination = as.character(doseCombination)) %>%
    left_join(., dose_probs_df, by = "doseCombination")
  
  # Switch between modes based on the provided flag
  if (mode == "interval") {
    # Interval mode: Categorise DLT probabilities into intervals
    patientData <- patientData %>%
      mutate(Interval = cut(DLT_prob, breaks = intervals, labels = interval_labels, right = TRUE))
  } else if (mode == "MTD") {
    # MTD mode: Determine if each dose combination is MTD
    patientData <- patientData %>%
      mutate(Interval = ifelse(DLT_prob >= lower_bound & DLT_prob <= upper_bound, "Yes", "No"))
  } else {
    stop("Invalid mode selected. Choose either 'interval' or 'MTD'.")
  }
  
  return(patientData)
}

average_summary_across_simulations <- function(simulation_results, prob_true_list, target_DLT_prob, indifference_interval, mode = 'MTD') {
  # Initialize a list to store summary results for each simulation
  summary_list <- list()
  
  # Process each simulation result
  for (i in seq_along(simulation_results)) {
    patientData <- simulation_results[[i]]$final_patient_data
    updated_patientData <- create_summary_table(prob_true_list, target_DLT_prob, indifference_interval, patientData, mode = mode)
    
    summary <- updated_patientData %>%
      group_by(Interval) %>%
      summarise(Number_of_Patients = n(), Number_of_DLTs = sum(outcome))
    
    summary_list[[i]] <- summary
  }
  # Combine all summaries into a single dataframe
  combined_summary <- bind_rows(summary_list, .id = "Simulation") %>%
    group_by(Interval) %>%
    summarise(Average_Number_of_Patients = mean(Number_of_Patients), Average_Number_of_DLTs = mean(Number_of_DLTs)) %>%
    ungroup()

  # Return the combined summary
  return(combined_summary)
}


average_summary <- average_summary_across_simulations(simulation_results, prob_true_list, 0.3, 0.1)
