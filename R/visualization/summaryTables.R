#' Create intervals around target DLT probability
#'
#' Creates intervals around the target DLT probability for categorizing dose combinations.
#'
#' @param target_DLT_prob Target DLT probability
#' @param interval_width Width of the interval
#'
#' @return Vector of interval boundaries
#' @export
create_intervals <- function(target_DLT_prob, interval_width) {
  lower_bound <- max(0, target_DLT_prob - interval_width / 2)
  upper_bound <- min(1, target_DLT_prob + interval_width / 2)
  
  # Define the intervals
  intervals <- c(0, lower_bound - interval_width, lower_bound, upper_bound, upper_bound + interval_width, 1)
  intervals <- unique(pmin(pmax(intervals, 0), 1))
  
  return(intervals)
}

#' Create summary table with categorization
#'
#' Creates a summary table of patient data with categorization by intervals or MTD status.
#'
#' @param prob_true_list List of true DLT probabilities
#' @param target_DLT_prob Target DLT probability
#' @param indifference_interval Width of indifference interval
#' @param patientData Patient data frame
#' @param mode Categorization mode ('interval' or 'MTD')
#'
#' @return Data frame with categorized patient data
#' @export
create_summary_table <- function(prob_true_list, target_DLT_prob, indifference_interval, patientData, mode = "MTD") {
  # Validate input
  if (!mode %in% c("interval", "MTD")) {
    stop("Invalid mode selected. Choose either 'interval' or 'MTD'.")
  }
  
  # Use the provided create_intervals function to generate intervals
  intervals <- create_intervals(target_DLT_prob, indifference_interval)
  interval_labels <- sapply(seq_along(intervals[-1]), function(i) {
    paste0("(", intervals[i], ",", intervals[i+1], "]")
  })
  
  # Define the lower and upper bounds for the MTD based on the interval around the target
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
    dplyr::mutate(doseCombination = as.character(doseCombination)) %>%
    dplyr::left_join(dose_probs_df, by = "doseCombination")
  
  # Switch between modes based on the provided flag
  if (mode == "interval") {
    # Interval mode: Categorise DLT probabilities into intervals
    patientData <- patientData %>%
      dplyr::mutate(
        Interval = cut(
          DLT_prob, 
          breaks = intervals, 
          labels = interval_labels, 
          right = TRUE
        )
      )
  } else if (mode == "MTD") {
    # MTD mode: Determine if each dose combination is MTD
    patientData <- patientData %>%
      dplyr::mutate(
        Interval = ifelse(
          DLT_prob >= lower_bound & DLT_prob <= upper_bound, 
          "Yes", 
          "No"
        )
      )
  }
  
  return(patientData)
}

#' Check if MTD falls within target interval
#'
#' Checks if the identified MTD falls within the acceptable interval around the target DLT probability.
#'
#' @param results Simulation results
#' @param prob_true_list List of true DLT probabilities
#' @param drugcombi_new Drug combination object
#' @param target_DLT_prob Target DLT probability
#' @param interval_width Width of indifference interval
#' @param quiet Suppress messages
#'
#' @return Logical indicating if MTD is within acceptable interval
#' @export
check_MTD_in_interval <- function(results, prob_true_list, drugcombi_new, target_DLT_prob, interval_width, quiet = TRUE) {
  # Retrieve the numeric index of the MTD from the results
  MTD_index <- results$RP2D
  
  # If MTD_index is NA, return NA
  if (is.na(MTD_index)) {
    if (!quiet) {
      message("MTD index is NA.")
    }
    return(NA)
  }
  
  # Convert the numeric index into dose label
  MTD_dose_label <- names(drugcombi_new$getDoseCombinationsLevel(MTD_index))
  
  # Extract the true probability of DLT for the identified MTD
  MTD_true_prob <- prob_true_list[[MTD_dose_label]]
  
  # Create intervals around the target DLT probability
  intervals <- create_intervals(target_DLT_prob, interval_width)
  
  # Check if the true probability of the MTD falls within the intervals
  within_interval <- MTD_true_prob >= intervals[3] && MTD_true_prob <= intervals[4]
  
  # If not quiet, print the message
  if (!quiet) {
    if (within_interval) {
      message(paste("The MTD falls within the acceptable interval:", intervals[3], "-", intervals[4]))
    } else {
      message(paste("The MTD does not fall within the acceptable interval:", intervals[3], "-", intervals[4]))
    }
  }
  
  # Return TRUE or FALSE
  return(within_interval)
}

#' Tabulate MTD selection frequency
#'
#' Calculates the percentage of simulations where each dose combination is identified 
#' as the Maximum Tolerated Dose (MTD) or Recommended Phase 2 Dose (RP2D).
#'
#' @param simulation_results A list of simulation result objects, each containing
#'   either an 'RP2D' element (single value) or an 'MTD' element (possibly multiple values)
#' @param nrow Integer specifying the number of rows in the dose matrix (typically the number of levels for drug A)
#' @param ncol Integer specifying the number of columns in the dose matrix (typically the number of levels for drug B)
#' @param type Character string indicating the type of analysis to perform:
#'   * "RP2D": Analyze the Recommended Phase 2 Dose selections (default)
#'   * "MTD": Analyze the Maximum Tolerated Dose selections (may include multiple values per simulation)
#'
#' @return A matrix of selection percentages, with dimensions nrow × ncol, 
#'   where each cell represents the percentage of simulations that identified 
#'   the corresponding dose combination as the MTD or RP2D
#'
#' @examples
#' \dontrun{
#' # Assuming simulation_results is a list of simulation outputs
#' # For a 3x4 dose combination matrix (drug A has 3 levels, drug B has 4 levels)
#' selection_matrix <- tabulate_MTD(simulation_results, nrow = 3, ncol = 4)
#' 
#' # For MTD analysis instead of RP2D
#' mtd_matrix <- tabulate_MTD(simulation_results, nrow = 3, ncol = 4, type = "MTD")
#' }
#'
#' @export
tabulate_MTD <- function(simulation_results, nrow, ncol, type = "RP2D") {
    # Input validation
    if (!is.list(simulation_results)) {
        stop("'simulation_results' must be a list")
    }
    if (!is.numeric(nrow) || !is.numeric(ncol) || nrow < 1 || ncol < 1) {
        stop("'nrow' and 'ncol' must be positive integers")
    }
    if (!type %in% c("RP2D", "MTD")) {
        stop("Invalid type. Use 'RP2D' or 'MTD'")
    }
    
    # Create a vector to store counts
    total_indices <- nrow * ncol
    counts <- rep(0, total_indices)
   
    # Length of the simulation list (denominator for percentages)
    total_simulations <- length(simulation_results)
    
    if (total_simulations == 0) {
        warning("Empty simulation results list provided")
        return(matrix(0, nrow = nrow, ncol = ncol))
    }
   
    if (type == "RP2D") {
        # Original RP2D counting logic
        for (sim in simulation_results) {
            rp2d <- sim$RP2D
            if (!is.na(rp2d) && rp2d >= 1 && rp2d <= total_indices) {
                counts[rp2d] <- counts[rp2d] + 1
            }
        }
       
        # Simple percentage for RP2D
        percentages <- (counts / total_simulations) * 100
       
    } else if (type == "MTD") {
        # MTD counting logic - handling multiple MTDs
        for (sim in simulation_results) {
            mtds <- sim$MTD
           
            if (!is.null(mtds) && length(mtds) > 0) {
                # Remove any NA or out of bounds values
                valid_mtds <- mtds[!is.na(mtds) & mtds >= 1 & mtds <= total_indices]
               
                if (length(valid_mtds) > 0) {
                    # Option 1: Simple percentage (counting each appearance)
                    counts[valid_mtds] <- counts[valid_mtds] + 1
                   
                    # Option 2: Weighted percentage (uncommment to use)
                    # weight <- 1/length(valid_mtds)
                    # counts[valid_mtds] <- counts[valid_mtds] + weight
                }
            }
        }
       
        # Calculate percentages based on total simulations
        percentages <- (counts / total_simulations) * 100
    }
   
    # Reshape into matrix
    result_matrix <- matrix(percentages, ncol = ncol, byrow = FALSE)
   
    return(result_matrix)
}