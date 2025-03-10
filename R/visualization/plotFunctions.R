#' Calculate cumulative DLT proportion
#'
#' Calculates cumulative dose-limiting toxicity proportion for each dose combination
#' up to a specified cohort.
#'
#' @param data Data frame containing patient outcomes
#' @param drugA Drug A object
#' @param drugB Drug B object
#' @param cohort_limit Maximum cohort to include
#' @param expand_grid Whether to include all possible combinations
#'
#' @return Data frame with cumulative DLT proportions
#' @export
calculate_cumulative_dlt <- function(data, drugA, drugB, cohort_limit, expand_grid = FALSE) {
  # Input validation
  if (!inherits(drugA, "Drug") || !inherits(drugB, "Drug")) {
    stop("drugA and drugB must be Drug objects")
  }
  
  # Filter data up to the specified cohort
  filtered_data <- data %>% 
    dplyr::filter(cohort <= cohort_limit)
  
  # Calculate cumulative DLT proportion for each dose combination level
  cumulative_dlt_data <- filtered_data %>%
    dplyr::arrange(cohort) %>%
    dplyr::group_by(doseCombination) %>%
    dplyr::summarise(
      cumulative_count = dplyr::n(),
      cumulative_DLT = sum(outcome),
      cumulative_DLT_proportion = cumulative_DLT / cumulative_count
    ) %>%
    dplyr::ungroup() 
  
  # Convert doseCombination to separate dose levels
  cumulative_dlt_data <- cumulative_dlt_data %>%
    tidyr::separate(doseCombination, into = c("Drug1", "Drug2"), sep = "\\.") %>%
    dplyr::mutate(
      Drug1 = factor(Drug1, levels = names(drugA$getDoseLevels())),
      Drug2 = factor(Drug2, levels = names(drugB$getDoseLevels()))
    )
  
  if (expand_grid) {
    all_combinations <- expand.grid(
      Drug1 = names(drugA$getDoseLevels()), 
      Drug2 = names(drugB$getDoseLevels())
    )
    
    cumulative_dlt_data <- all_combinations %>% 
      dplyr::left_join(cumulative_dlt_data, by = c("Drug1", "Drug2"))
  }
  
  return(cumulative_dlt_data)
}

#' Prepare step data for visualization
#'
#' Prepares data for plotting step functions in the dose allocation visualization.
#'
#' @param sim_output Simulation output
#' @param drugA Drug A object
#' @param drugB Drug B object
#' @param cohort_limit Maximum cohort to include
#'
#' @return Data frame with step function coordinates
#' @export
prepare_step_data <- function(sim_output, drugA, drugB, cohort_limit) {
  # Add toxicity status to the dataframe
  pipe_estimate <- sim_output$simulation_results[[cohort_limit]]$best_config
  nlevel_drugA <- length(drugA$getDoseLevels())
  nlevel_drugB <- length(drugB$getDoseLevels())

  xlevel <- 1:(nlevel_drugA + 1) - 0.5
  ylevel <- c(
    apply(
      matrix(pipe_estimate, ncol = nlevel_drugA), 
      2, 
      function(i) {
        min(which(i == 0), nlevel_drugB + 1) - 0.5
      }
    ), 
    0.5
  )

  # Prepare the step data
  step_data <- data.frame(
    x = xlevel,
    y = ylevel
  )

  return(step_data)
}

#' Plot DLT proportion
#'
#' Creates a visualization of DLT proportions across dose combinations.
#'
#' @param sim_output Simulation output
#' @param drugA Drug A object
#' @param drugB Drug B object
#' @param cohort_seq Cohort sequence to visualize
#'
#' @return ggplot object
#' @import ggplot2
#' @export
plot_dlt_proportion <- function(sim_output, drugA, drugB, cohort_seq) {
  # Calculate cumulative DLT and prepare step data
  cumulative_dlt_df <- calculate_cumulative_dlt(
    sim_output$final_patient_data, 
    drugA, 
    drugB, 
    cohort_seq, 
    expand_grid = TRUE
  )
  
  step_data <- prepare_step_data(sim_output, drugA, drugB, cohort_seq)
  n_doses_drugA <- length(drugA$getDoseLevels())
  n_doses_drugB <- length(drugB$getDoseLevels())
  
  # Create the plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = cumulative_dlt_df, 
      ggplot2::aes(x = Drug1, y = Drug2, fill = cumulative_DLT_proportion), 
      alpha = 0.5, 
      size = 2
    ) +
    ggplot2::scale_x_discrete(
      breaks = names(drugA$getDoseLevels()), 
      expand = c(0, 0), 
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(
      breaks = names(drugB$getDoseLevels()), 
      expand = c(0, 0), 
      drop = FALSE
    ) +
    ggplot2::geom_step(
      ggplot2::aes(x = x, y = y), 
      data = step_data, 
      linewidth = 1.5, 
      color = "black", 
      linetype = 2
    ) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
      limits = c(0, 1),
      na.value = 'white',
      name = "DLT Proportion"
    ) +
    ggplot2::geom_text(
      data = cumulative_dlt_df %>% dplyr::filter(!is.na(cumulative_DLT_proportion)), 
      ggplot2::aes(
        x = as.numeric(Drug1), 
        y = as.numeric(Drug2), 
        label = paste(cumulative_count, "(", cumulative_DLT, ")", sep = "")
      ), 
      color = "black"
    ) +
    ggplot2::labs(x = "Drug 1", y = "Drug 2") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}
#' Plot PIPE adaptation across cohorts
#'
#' Creates a grid of plots showing how the PIPE algorithm adapts its dose recommendations
#' as new patient outcomes are observed across multiple cohorts in a single simulation.
#' Each plot displays the cumulative DLT proportions and the current boundary
#' between estimated safe and unsafe doses after each cohort.
#'
#' @param simulation_result A single simulation result object
#' @param drugA Drug A object
#' @param drugB Drug B object
#' @param max_cohorts Maximum number of cohorts to display, or NULL to use all available cohorts
#' @param ncol Number of columns in the grid
#' @param title_text Title text for the plot
#' @param save_path File path to save the plot, or NULL to not save
#' @param width Width of the saved plot in inches
#' @param height Height of the saved plot in inches
#'
#' @return A ggplot object with the combined plots
#' @import ggplot2
#' @import cowplot
#' @export

plot_pipe_adaptation <- function(simulation_result, drugA, drugB, max_cohorts = NULL,
                                     ncol = 4, title_text = "Cumulative DLT Proportion by Dose Combination",
                                     save_path = NULL, width = 30, height = 40) {
  
  # Input validation
  if (is.null(simulation_result)) {
    stop("Simulation result must be provided")
  }
  if (is.null(drugA) || is.null(drugB)) {
    stop("Drug objects must be provided")
  }
  
  # Determine the maximum cohort number from simulation_result
  if (is.null(max_cohorts)) {
    if (!is.null(simulation_result$final_patient_data)) {
      max_cohorts <- max(simulation_result$final_patient_data$cohort, na.rm = TRUE)
    } else if (!is.null(simulation_result$simulation_results)) {
      max_cohorts <- length(simulation_result$simulation_results)
    } else {
      stop("Cannot determine maximum cohort number from simulation_result")
    }
  }
  
  # Create sequence of cohorts to plot
  cohort_seqs <- 1:max_cohorts
  
  # Create list of plots for each cohort
  plots <- lapply(cohort_seqs, function(cohort_seq) {
    plot_dlt_proportion(simulation_result, drugA, drugB, cohort_seq)
  })
  
  # Extract the legend from the first plot
  legend <- cowplot::get_legend(plots[[1]] + 
                                theme(legend.position = "bottom",
                                      legend.box.margin = margin(0, 0, 0, 0)))
  
  # Remove legends from all plots
  plots <- lapply(plots, function(plot) {
    plot + theme(legend.position = "none")
  })
  
  # Combine plots into a grid using cowplot
  combined_plot <- cowplot::plot_grid(
    plotlist = plots, 
    ncol = ncol, 
    align = 'v'
  )
  
  # Add a common title
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(title_text, 
                        fontface = 'bold', 
                        size = 14)
  
  # Combine title, plots, and legend
  final_plot <- cowplot::plot_grid(
    title, 
    combined_plot, 
    legend, 
    ncol = 1, 
    rel_heights = c(0.1, 1, 0.1)
  )
  
  # Save the plot if a path is provided
  if (!is.null(save_path)) {
    cowplot::save_plot(
      filename = save_path,
      plot = final_plot,
      base_height = height,
      base_width = width
    )
  }
  
  return(final_plot)
}

#' Animate PIPE adaptation process
#'
#' Creates an animation showing how the PIPE algorithm dynamically adapts its
#' dose recommendations as new patient outcomes are observed throughout the trial.
#' The animation shows the progression of DLT proportions and the evolving boundary
#' between safe and unsafe doses.
#' 
#' @param simulation_result A single simulation result object
#' @param drugA Drug A object
#' @param drugB Drug B object
#' @param max_cohorts Maximum number of cohorts to include in the animation, or NULL to use all available
#' @param transition_length Length of transitions between states
#' @param state_length Length of time to display each state
#' @param fps Frames per second in the animation
#' @param save_path File path to save the animation, or NULL to not save
#' @param width Width of the animation in pixels
#' @param height Height of the animation in pixels
#'
#' @return A gganim object
#' @import ggplot2
#' @import gganimate
#' @export

animate_pipe_adaptation <- function(simulation_result, drugA, drugB, max_cohorts = NULL,
                                transition_length = 2, state_length = 1, fps = 10,
                                save_path = "dlt_proportion_animation.gif",
                                width = 800, height = 600) {
  
  # Input validation
  if (is.null(simulation_result)) {
    stop("Simulation result must be provided")
  }
  if (is.null(drugA) || is.null(drugB)) {
    stop("Drug objects must be provided")
  }
  
  # Determine the maximum cohort number from simulation_result
  if (is.null(max_cohorts)) {
    if (!is.null(simulation_result$final_patient_data)) {
      max_cohorts <- max(simulation_result$final_patient_data$cohort, na.rm = TRUE)
    } else if (!is.null(simulation_result$simulation_results)) {
      max_cohorts <- length(simulation_result$simulation_results)
    } else {
      stop("Cannot determine maximum cohort number from simulation_result")
    }
  }
  
  # Create sequence of cohorts for the animation
  cohort_seqs <- 1:max_cohorts
  
  # Prepare data for animation
  all_data <- data.frame()
  all_steps <- data.frame()
  
  # Collect data for all cohorts
  for (cohort in cohort_seqs) {
    # Calculate cumulative DLT
    cohort_data <- calculate_cumulative_dlt(
      simulation_result$final_patient_data, 
      drugA, 
      drugB, 
      cohort, 
      expand_grid = TRUE
    )
    cohort_data$cohort <- cohort
    
    # Prepare step data
    step_data <- prepare_step_data(simulation_result, drugA, drugB, cohort)
    step_data$cohort <- cohort
    
    # Combine data
    all_data <- rbind(all_data, cohort_data)
    all_steps <- rbind(all_steps, step_data)
  }
  
  # Create base plot
  base_plot <- ggplot() +
    geom_tile(
      data = all_data, 
      aes(x = Drug1, y = Drug2, fill = cumulative_DLT_proportion), 
      alpha = 0.5, 
      size = 2
    ) +
    scale_x_discrete(
      breaks = names(drugA$getDoseLevels()), 
      expand = c(0, 0), 
      drop = FALSE
    ) +
    scale_y_discrete(
      breaks = names(drugB$getDoseLevels()), 
      expand = c(0, 0), 
      drop = FALSE
    ) +
    geom_step(
      aes(x = x, y = y, group = cohort), 
      data = all_steps, 
      linewidth = 1.5, 
      color = "black", 
      linetype = 2
    ) +
    scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
      limits = c(0, 1),
      na.value = 'white',
      name = "DLT Proportion"
    ) +
    geom_text(
      data = all_data %>% filter(!is.na(cumulative_DLT_proportion)), 
      aes(
        x = as.numeric(Drug1), 
        y = as.numeric(Drug2), 
        label = paste(cumulative_count, "(", cumulative_DLT, ")", sep = "")
      ), 
      color = "black"
    ) +
    labs(
      x = "Drug 1", 
      y = "Drug 2",
      title = "Cohort: {closest_state}"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Create animation
  anim <- base_plot + 
    gganimate::transition_states(
      cohort,
      transition_length = transition_length,
      state_length = state_length
    ) +
    gganimate::enter_fade() + 
    gganimate::exit_fade()
  
  # Save the animation if a path is provided
  if (!is.null(save_path)) {
    gganimate::anim_save(
      filename = save_path,
      animation = anim,
      width = width,
      height = height,
      fps = fps
    )
  }
  
  return(anim)
}