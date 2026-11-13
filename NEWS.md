# gPIPE (development version)

* Requires R >= 4.4.
* CI: R CMD check on R-devel (non-blocking, weekly), release and oldrel-1 for `main` and `dev`; Codecov coverage reporting; lintr on changed files (#16).
* New S7 class `DoseGrid` for dose combinations of any number of drugs, with helpers to convert between the array and vector views (#18).
* New S7 API (#18): `PipeDesign`, with the `ms_design()` (Mander & Sweeting 2015) and `cheung_design()` (Cheung et al. 2022) constructors; `estimate()`, `next_dose()` and a `simulate()` method. One exact estimator (a minimum cut) is used for any number of drugs. The Mander & Sweeting safety probabilities are computed exactly without listing the monotone contours, so the safety rule also works for three drugs (#31).
* Faster simulation: the trial loop reads the grid from a plain list made once per design, and no longer runs `match.arg()` at every step. Results are unchanged (#31).
* The minimum cut of the estimator now runs in C++ through Rcpp (a new dependency), so the package needs a compiler to install from source. A 5 x 5 trial of 125 patients simulates about four times faster, and grids of thousands of combinations become practical. Results are unchanged (#33).
* New `selection` options for the Cheung rule: `"fewest"` (equal probability on the candidates with the fewest patients), `"inverse_n"` (proportional to 1 / (n + 1)) and `"interval"` (proportional to the posterior probability that the parameter is within `delta` of `theta`; binary outcome only). The sampling set, the defaults and the results of `"equal"` and `"prob"` are unchanged (#45).
* New `shift` property of `PipeDesign` for `rule = "cheung"`: decisions use the threshold `theta + shift`, and `theta` stays the clinical target (Riviere's pick, `summary()` and `calibrate_eps()` still use `theta`). With `eps` < 0.5 a dose at exactly `theta` is called safe in only about a fraction `eps` of trials, and a positive `shift` corrects this. The default 0 changes nothing (#45).
* New option `startup = "diagonal"` for `PipeDesign` (Cheung rule, binary outcome, MTD): the diagonal start-up phase of Riviere et al. (2014, Section 2.1.3). The trial raises every drug that is not at its top level by one level per cohort until the first toxicity, a top combination or the edge of a subset grid, then switches to the usual rule. `next_dose()` follows the phase, and the simulation summary has a new column `n_startup`. The default `"none"` is unchanged (#51).
* Fixed: per-dose priors no longer give `NA` posteriors; trials no longer stop silently when no border dose is in the neighbourhood; the taper reaches `eps` at the last patient.
* Fixed: the recommended phase II doses of a Mander & Sweeting design are now the closest doses among those that pass the safety rule, as in pipe.design. Before, a dose that was safe on the contour but failed the safety rule stopped the dose below it from being recommended (about 1 trial in 450 in the Table III settings) (#31).
* New `pick = "isotonic"` option of `PipeDesign` for the single recommended dose (binary outcome, MTD): the tried dose whose weighted isotonic estimate of the toxicity rate is closest to the target, computed with the min-cut solver. The default `"interval"` is Riviere's rule as before (#54).
* New `selection = "closest"` for the Cheung rule (binary outcome, MTD): the next cohort goes to the dose in the sampling set whose isotonic estimate of the toxicity rate is closest to `theta` (an untried dose counts as `theta`; ties share equal probability). The other selections and the defaults are unchanged (#57).
* Golden tests against pipe.design 0.5.1 decision snapshots and Cheung et al. (2022) Table 1, run per implementation so each refactor step enables its own tests (#17).

## Breaking changes

* The reference-class API is removed, together with its exports:
  * classes `Drug`, `DrugCombination`, `DoseConfiguration`, `DrugCombinationModel`, `PipeEstimator`, `PatientDataModel`, `AdmissibleCombinationRule`, `AdjacentAdmissible`, `BorderAdmissible`, `ClosestAdmissible`, `NeighbourAdmissible`, `SafetyAdmissible`, `SelectionStrategy`, `EqualRandomisationStrategy`, `InverseDistanceStrategy`, `InverseSqrtSampleSizeStrategy`, `PosteriorProbabilityStrategy`, `SmallestSampleSizeStrategy`;
  * functions `createDrug()`, `createDrugCombi()`, `createDoseConfig()`, `createDrugCombinationModel()`, `createPipeEstimator()`, `createPatientDataModel()`, `runTrialSimulation()`, `calibrate_epsilon()`, `monotonic_matrices()`, `generate_gamma_config()`, `calculateNeighbourSum()`, `calculate_cumulative_dlt()`, `plot_dlt_proportion()`, `plot_pipe_adaptation()`, `prepare_step_data()`, `create_intervals()`, `create_summary_table()`, `tabulate_MTD()`, `check_MTD_in_interval()`.
* Calibration, plots and summary tables come back on the S7 API (#29).
* Dependencies dplyr, tidyr, ggplot2, cowplot, RColorBrewer, methods and utils are dropped.

# gPIPE 0.1.0

## gPIPE 0.1.0

### New Features

* Initial implementation of the generalized PIPE methodology with epsilon-tapering
* Comprehensive class structure for dose-finding study components:
  * `DrugCombination` class for managing drug dose combinations
  * `DoseConfiguration` class for handling dose configurations
  * `PatientDataModel` class for patient data and trial progression
  * `PipeEstimator` class implementing the core PIPE methodology
* Multiple admissibility rules:
  * Safety-based rule
  * Neighbor-based rule
  * Border-based rule
  * Closest boundary rule
  * Adjacent dose rule
* Various dose selection strategies:
  * Smallest sample size strategy
  * Posterior probability strategy
  * Inverse square root sample size strategy
  * Inverse distance strategy
  * Equal randomization strategy
* Utility functions for:
  * Generating monotonic matrices
  * Calculating neighbor sums
  * Epsilon calibration
  * Finding boundary regions
* Comprehensive documentation and examples

### Internal Changes

* Established reference class inheritance hierarchies for extensibility
* Implemented modular structure for combining different admissibility rules and selection strategies
* Added roxygen2 documentation for all classes and functions
