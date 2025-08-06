#' Generate Valid Gamma Configurations
#'
#' Generates all valid gamma configurations (binary vectors) that satisfy
#' monotonicity constraints for dose-finding studies.
#'
#' @param n_dose Integer indicating the number of dose levels
#'
#' @return A data frame containing all valid gamma configurations where
#'   each row represents a valid configuration and each column represents
#'   a dose level (1 = acceptable, 0 = unacceptable)
#'
#' @details
#' This function generates all possible binary configurations for n_dose levels
#' and filters them to keep only those that satisfy the monotonicity constraint
#' (non-decreasing sequences). This ensures that if a dose level is acceptable,
#' all lower dose levels are also acceptable.
#'
#' @examples
#' # Generate valid configurations for 3 dose levels
#' configs <- generate_gamma_config(3)
#' print(configs)
#'
#' @importFrom dplyr slice
#' @export
generate_gamma_config <- function(n_dose) {
    # Generate all possible binary combinations
    gamma_grid <- expand.grid(rep(list(0:1), n_dose))
    
    # Find indices of configurations that satisfy monotonicity (non-decreasing)
    valid_indices <- which(apply(gamma_grid, 1, function(x) all(diff(x) >= 0)))
    
    # Select only the valid configurations
    gamma_grid_valid <- dplyr::slice(gamma_grid, valid_indices)
    
    return(gamma_grid_valid)
}
