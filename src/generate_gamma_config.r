generate_gamma_config <- function(n_dose) {
    gamma_grid <- rep(list(0:1), n_dose) %>%
        expand.grid() %>%
        as.data.frame()
    gamma_grid_valid <- gamma_grid %>%
        with(which(apply(., 1, function(x) all(diff(x) >= 0)))) %>%
        slice(gamma_grid, .)

    return(gamma_grid_valid)
}
