# Objective function to minimize
# Here, params is a vector of two elements: params[1] = a, params[2] = b
objective_function <- function(params, target_median, prior_ss) {
  a <- params
  median_estimate <- qbeta(0.5, a, prior_ss - a)
  return(median_estimate - target_median)
}

# Main function to find a and b
find_beta_parameters <- function(prior_med, prior_ss) {
  results <- lapply(1:length(prior_med), function(i) {
    target_median <- prior_med[i]
    # Initial guesses for a and b based on median and sample size
    initial_guess <- prior_ss[i] * target_median
    
    # Optimizing
    optim_res <- uniroot(objective_function, target_median = target_median, prior_ss = prior_ss[i], c(1e-4, 1))
    
    return(c(a = optim_res$root, b = prior_ss[i] - optim_res$root))
  })
  
  # Convert results to a suitable format
  a_values <- sapply(results, function(x) x["a"])
  b_values <- sapply(results, function(x) x["b"])
  
  return(list(a = a_values, b = b_values))
}
