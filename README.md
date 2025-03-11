
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gPIPE

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/gPIPE)](https://CRAN.R-project.org/package=gPIPE)
[![R-CMD-check](https://github.com/xlai/gPIPE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xlai/gPIPE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Generalized PIPE Design with Epsilon-Tapering for Dose-Finding Studies

The gPIPE package implements a generalized Product of Independent Beta
Probabilities Escalation (PIPE) design incorporating epsilon-tapering
for dual-agent dose-finding studies. Building on Cheung et al.’s (2022)
asymmetric posterior gain and epsilon-tapering framework, gPIPE improves
dose-space exploration while effectively controlling false discovery
rates. The implementation features flexible admissibility rules,
multiple dose selection strategies, and robust simulation tools for
design calibration, evaluation, and comparison with existing
dose-finding methodologies.

## Installation

You can install the development version of gPIPE from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("xlai/gPIPE")
```

## Example

Here’s a basic example of setting up a dose-finding trial using gPIPE:

``` r
library(gPIPE)

# Create a drug combination trial
trial <- DrugCombiConstructor$new(
  drug_A = list(name = "Drug A", dose_levels = c(1, 3, 5, 10, 15)),
  drug_B = list(name = "Drug B", dose_levels = c(0.1, 0.5, 1.5, 3, 5))
)

# Setup PIPE model with epsilon tapering
pipe_model <- PIPEConstructor$new(
  trial = trial,
  target_dlt_rate = 0.3,
  epsilon = 0.05
)

# Run a simulation
results <- runTrialSimulation(
  pipe_model = pipe_model,
  n_simulations = 100,
  cohort_size = 3,
  max_sample_size = 36
)
```

## Features

- Implementation of the generalized PIPE design for dose-finding studies
- Support for epsilon-tapering to improve dose-space exploration
- Flexible admissibility rules for dose combinations
- Multiple dose selection strategies
- Simulation tools for design calibration and evaluation

## Documentation

For more detailed information, please refer to the package vignette:

``` r
vignette("gPIPE")
```

## License

This package is licensed under the
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html).
