# gPIPE (development version)

* Requires R >= 4.4.
* CI: R CMD check on R-devel (non-blocking, weekly), release and oldrel-1 for `main` and `dev`; Codecov coverage reporting; lintr on changed files (#16).
* New S7 class `DoseGrid` for dose combinations of any number of drugs, with helpers to convert between the array and vector views (#18).
* New S7 API (#18): `PipeDesign`, with the `ms_design()` (Mander & Sweeting 2015) and `cheung_design()` (Cheung et al. 2022) constructors; `estimate()`, `next_dose()` and a `simulate()` method. Exact solvers (enumeration and min-cut) work on any number of drugs.
* Fixed: per-dose priors no longer give `NA` posteriors; trials no longer stop silently when no border dose is in the neighbourhood; the taper reaches `eps` at the last patient.
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
