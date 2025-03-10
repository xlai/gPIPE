map_vector <- function(input_vector, target_dlt_rate, interval_width) {
  # Get equally spaced out probabilities set by target_dlt_rate and interval_width
  p_prior <- dfcrm::getprior(interval_width / 2, target_dlt_rate, 2, 3)
  
  # Define the lower and upper bounds of the target interval
  lower_bound <- target_dlt_rate - interval_width / 2
  upper_bound <- target_dlt_rate + interval_width / 2
  
  # Initialize the output vector with NA (or any placeholder)
  output_vector <- numeric(length(input_vector))
  
  # Apply the mapping rules explicitly
  output_vector[input_vector < lower_bound] <- p_prior[1]
  output_vector[input_vector >= lower_bound & input_vector <= upper_bound] <- p_prior[2]
  output_vector[input_vector > upper_bound] <- p_prior[3]
  
  return(output_vector)
}

map_vector2 <- function(input_vector, target_dlt_rate, odds_ratio = 2) {
  # Define the lower and upper bounds of the target interval
  lower_bound <- target_dlt_rate / (odds_ratio + target_dlt_rate*( 1- odds_ratio))
  upper_bound <- target_dlt_rate * odds_ratio / ( 1 - target_dlt_rate*( 1- odds_ratio))
  
  # Initialize the output vector with NA (or any placeholder)
  output_vector <- numeric(length(input_vector))
  
  # Apply the mapping rules explicitly
  output_vector[input_vector < target_dlt_rate] <- lower_bound
  output_vector[input_vector >= target_dlt_rate] <- upper_bound
  
  return(output_vector)
}