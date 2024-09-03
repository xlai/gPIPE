
epsilon_tapering <- function(n_total, n_current, epsilon, taper_type = "linear", custom_taper = NULL) {
  
  if (!is.null(custom_taper)) {
    # Use the custom tapering function provided by the user
    epsilon_current <- custom_taper(n_total, n_current, epsilon)
  } else {
    if (taper_type == "linear") {
      # Linear tapering
      epsilon_current <- (n_total - n_current) / (n_total + 1) * 0.5 + 
                         (n_current + 1) / (n_total + 1) * epsilon
    } else if (taper_type == "quadratic") {
      # Quadratic tapering: slow down tapering at the beginning
      alpha <- n_current / n_total
      epsilon_current <- (1 - alpha)^2 * 0.5 + alpha^2 * epsilon
    } else {
      stop("Unsupported taper_type. Please use 'linear', 'quadratic', or provide a custom taper function.")
    }
  }
  
  return(epsilon_current)
}