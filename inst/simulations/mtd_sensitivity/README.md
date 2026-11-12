# MTD component sensitivity study (issue #45, steps 1 and 2)

How much does each component of PIPE with epsilon-tapering matter for
finding the maximum tolerated dose (MTD) contour with binary toxicity? This
folder holds the scaffold and the first runs: the baseline and two reference
lines. Later steps add configurations that change one component at a time.

- `common.R`: the scenario, the configuration table `mtd_configs()`, the
  design built from a row, and the summary of one simulation.
- `calibrate.R`: calibrates epsilon for the baseline design.
- `run.R`: simulates configurations and stores the summaries.
- `report.R`: prints the stored results.
- `results/`: the stored runs (see below).

To add a configuration, add one row in `mtd_configs()` (with a new `id` and
a new `seed`) and run `Rscript run.R 4 <id>`.

## Scenario

Riviere et al. (2014), scenario 1: drug A with 5 levels, drug B with 3
levels, target theta = 0.30. True toxicity (rows are drug A levels):

| | B1 | B2 | B3 |
|---|---|---|---|
| A1 | 0.05 | 0.10 | 0.15 |
| A2 | 0.10 | 0.15 | 0.30 |
| A3 | 0.15 | 0.30 | 0.45 |
| A4 | 0.30 | 0.45 | 0.50 |
| A5 | 0.45 | 0.55 | 0.60 |

The true MTDs are the three cells with toxicity 0.30: A4.B1, A3.B2, A2.B3.
The values and the combination order (A varies fastest) were checked
against the scenario file and the dose grid.

## Baseline

`cheung_design()`, target MTD, binary outcome, epsilon-tapering
(`taper = TRUE`), weights `n+1`, selection `equal`, N = 60 in cohorts of 3,
start at combination 1, single pick with delta = 0.1 inside the final set.
The prior is flat: Beta with median 0.30 at every combination and strength
a + b = 1/15 (1/K). The median is solved exactly (`method = "exact"`).

Epsilon is calibrated by `calibrate_eps()` with alpha = 0.05 and band = 0.1:
all 56 monotone configurations of the 5 x 3 grid (toxic doses at theta + 0.1,
safe doses at theta - 0.1), 600 trials per configuration and evaluation,
bisection on [0.01, 0.5] down to a bracket of 0.005, seed 1. The path of
every epsilon tried (FDR, TPR, SE) is in `results/calibration_baseline.csv`,
with the settings. `run.R` reads epsilon from that file.

## Other configurations

- `eps_0.5`: the baseline with epsilon = 0.5 (tapering has no effect).
- `ref_ms`: the Mander & Sweeting design, `ms_design()` with its defaults
  (closest, neighbouring, safety 0.8, `min_S`), the same flat prior
  (median 0.30, strength 1/15, found by bisection as `ms_design()` does),
  N = 60, cohorts of 3.
- Riviere's published PCS for scenario 1 (`results/riviere2014_scenario1_pcs.csv`)
  is a reference line only. Those designs are different and they used 2000
  replications.

## Metrics

Each configuration runs 10000 trials. Every value has a Monte Carlo
standard error. Whole-trial metrics from `summary()` with band = 0.1:

- `fdr`, `tpr`: Cheung et al. (2022), eqs. 7 and 8. Proportions.
- `pcs`: proportion of trials whose single pick is a true MTD
  (toxicity = theta), as in Riviere et al. A trial with no pick counts as
  a failure. It is the sum of `pct_pick` over the three true MTDs, divided
  by 100.
- `pct_patients_above`, `pct_patients_in_band`, `pct_patients_below`:
  percentage of patients treated at combinations with toxicity at least
  0.40, between 0.20 and 0.40, and at most 0.20.
- `p_set_false`: proportion of trials whose final set has a combination
  with toxicity at least 0.40.
- `mean_recommended`: mean size of the final set.
- `pct_stopped`, `mean_patients`: percentage of trials stopped early, and
  mean patients enrolled.

Per combination (one row for each of the 15):

- `pct_selected`: percentage of trials with the combination in the final
  set.
- `pct_pick`: percentage of trials whose single pick is the combination.
- `mean_n`: mean patients treated.

Mechanism metrics, from the dose of each cohort:

- `first_contour_cohort`: number of the first cohort (the first is 1) given
  a combination on the true contour (toxicity at least theta). Trials that
  never do this are left out of this mean, and their share is
  `pct_never_contour`.
- `n_distinct`: number of different combinations given.
- `pct_skip`: percentage of trials with a skip. A skip is a cohort, after
  the first, whose combination is more than one level above every earlier
  combination in some drug. It is not a skip when the combination is
  within one level above, in every drug, of at least one earlier one. A
  lower level in a drug counts as within one level.
- `n_dlt`: mean number of DLTs per trial.
- `pct_set_has_mtd`: percentage of trials whose final set has at least one
  true MTD.
- `n_mtd_tried`: mean number of true MTDs given to at least one patient.

Units: `fdr`, `tpr`, `pcs`, `p_set_false` are proportions. Names that start
with `pct_` are percentages. The others are counts.

## Results files

- `mtd_sensitivity.csv`: long format: `id`, `metric`, `combo` (empty for a
  whole-trial metric), `estimate`, `mc_se`.
- `mtd_sensitivity_runs.csv`: per configuration, trials, seed, epsilon used,
  and seconds.
- `calibration_baseline.csv`, `riviere2014_scenario1_pcs.csv`.

## Seeds

Each configuration has its own `seed` in `mtd_configs()`. It runs its 10000
trials in one process, with `set.seed(seed)`. The parallel loop is over
configurations, so results do not depend on the number of cores or on which
other configurations run. The calibration uses seed 1.

## How to run

From this folder:

    Rscript calibrate.R          # about 10 minutes, one core
    Rscript run.R 4              # all configurations, about 4 minutes on 4 cores
    Rscript run.R 4 eps_0.5      # only these, and only their rows change
    Rscript report.R [id]        # tables; grids for one configuration

Run times for 10000 trials: baseline and `eps_0.5` about 26 seconds each,
`ref_ms` about 230 seconds (the safety rule is slower). `run.R` needs the
installed gPIPE.

`tests/testthat/test-simulation-mtd-sensitivity.R` checks the stored results.
