# gPIPE (development version)

* Requires R >= 4.4.
* CI: R CMD check on R-devel (non-blocking, weekly), release and oldrel-1 for `main` and `dev`; Codecov coverage reporting; lintr on changed files (#16).
* Golden tests against pipe.design 0.5.1 decision snapshots and Cheung et al. (2022) Table 1, run per implementation so each refactor step enables its own tests (#17).

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
