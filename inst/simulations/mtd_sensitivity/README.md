# MTD component sensitivity study (issue #45)

How much does each component of PIPE with epsilon-tapering matter for
finding the maximum tolerated dose (MTD) contour with binary toxicity? This
folder holds the baseline, two reference lines, rows that change one
factor at a time from the baseline, rows for the selection rules and the
band shift, and a few exploratory rows that change several things at once.

- `common.R`: the scenario, the configuration table `mtd_configs()`, the
  design built from a row, and the summary of one simulation.
- `calibrate.R`: calibrates epsilon for the baseline design, for a given
  alpha, band and shift.
- `run.R`: simulates configurations and stores the summaries.
- `report.R`: prints the stored results by factor, with the difference from
  the baseline, and a table of all rows sorted by PCS.
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

Two more calibrations use alpha = 0.10 and 0.20 and nothing else changed
(`results/calibration_alpha_0.10.csv`, `results/calibration_alpha_0.20.csv`).
A configuration names the file its epsilon comes from in its `calibration`
column.

## Other configurations

- `eps_0.5`: the baseline with epsilon = 0.5 (tapering has no effect).
- `ref_ms`: the Mander & Sweeting design, `ms_design()` with its defaults
  (closest, neighbouring, safety 0.8, `min_S`), the same flat prior
  (median 0.30, strength 1/15, found by bisection as `ms_design()` does),
  N = 60, cohorts of 3.
- Riviere's published PCS for scenario 1 (`results/riviere2014_scenario1_pcs.csv`)
  is a reference line only. Those designs are different and they used 2000
  replications.

## One-factor-at-a-time rows

Each row changes one thing from the baseline. **Every row keeps the baseline
epsilon (0.113) except the epsilon rows, the calibration rows and the band shift rows.** Epsilon
is not recalibrated when N, the prior, the cohort size or any other factor
changes. So a change in N or in the prior also changes the false discovery
rate that epsilon was chosen to control, and the rows are not at equal FDR.
Compare the FDR column as well as PCS.

| factor | rows |
|---|---|
| eps | `eps_0.4`, `eps_0.3`, `eps_0.15`, `eps_0.05` (and `eps_0.5`) |
| calibration alpha | `alpha_0.10` (epsilon 0.217), `alpha_0.20` (epsilon 0.492) |
| taper | `taper_off`: `taper = FALSE` |
| weights | `weights_1`: weights `"1"` |
| prior strength | `prior_strength_1`, `prior_strength_4`: flat, median 0.30, strength 1 and 4 |
| prior form | `prior_beta_1_1`: Beta(1,1). `prior_beta_0.5_0.5`: Beta(0.5,0.5). `prior_correct`: increasing medians equal to the truth of scenario 1. `prior_too_low`: medians equal to the truth of Riviere scenario 3 (the prior thinks doses are less toxic than they are). `prior_too_high`: medians equal to scenario 2 (more toxic). The three increasing priors have strength 1/15. |
| prior form x strength | `prior_correct_strength_1`: the correct medians with strength 1 |
| selection | `selection_prob`: selection `"prob"` |
| cohort | `cohort_1`: cohorts of 1 |
| N | `N_30`, `N_120` |
| start | `start_top`: start at A5.B3 (combination 15) |
| single pick | `single_pick_off`: `within_set = FALSE` |

## Selection, band shift and combination rows

**Selection.** `selection_fewest`, `selection_inverse_n` and
`selection_interval` use the other ways to draw the next dose from the
sampling set (`selection = "fewest"`, `"inverse_n"`, `"interval"`). They keep
the baseline epsilon.

**Band shift.** The shift moves the decision threshold. The design decides
with P(pi > theta + shift) instead of P(pi > theta). Theta stays the target:
the single pick, `summary()` and `calibrate_eps()` still use 0.30. With a
shift, a dose whose true toxicity is exactly 0.30 is called safe more often,
so true MTDs are kept in the final set more often.

The shift rows take the band b of the calibration as the shift. For
b = 0.1, 0.075, 0.05 and 0.025, the design has `shift = b`, and epsilon is
calibrated again for it with alpha = 0.05 and band = b over the same 56
configurations. The simulated truths sit at theta +- b, so the toxic ones are
exactly on the decision threshold theta + b. The prior stays the baseline's
(median 0.30, strength 1/15). The calibrations are in
`results/calibration_shift_<b>.csv` (columns `band` and `shift` record what
they used). `summary()` is still called with band = 0.1 for every row, so
FDR, TPR and PCS are on the same footing as the baseline.

| factor | rows |
|---|---|
| selection | `selection_fewest`, `selection_inverse_n`, `selection_interval` |
| band shift | `shift_0.1`, `shift_0.075`, `shift_0.05`, `shift_0.025`: shift b with epsilon calibrated at band b |
| band shift, fixed eps | `shift_0.1_eps_fixed`: shift 0.1 with the baseline epsilon 0.113, to separate the effect of the shift from the effect of the new epsilon |
| combination | `eps_0.4_pick_all`, `eps_0.3_pick_all` (epsilon 0.4 and 0.3 with `within_set = FALSE`); `shift_0.1_pick_all`, `shift_0.05_pick_all` (shift with its calibrated epsilon and `within_set = FALSE`); `shift_0.1_pick_all_fewest`, `shift_0.1_pick_all_interval` (as `shift_0.1_pick_all`, with selection `"fewest"` or `"interval"`) |

The combination rows are exploratory. Each changes more than one thing from
the baseline, so a difference cannot be put down to one factor. The
`shift` column of the configuration table is 0 for every other row.

The two other scenario truths are in `common.R` (`mtd_truth_s2`,
`mtd_truth_s3`), copied from `riviere2014_scenarios_all14.csv` in the same
combination order. A Beta prior with median 0.5 and strength 2 is Beta(1,1),
and strength 1 is Beta(0.5,0.5), so both are flat rows of the table. The
prior column `prior_form` is `flat` or the scenario whose truth gives the
medians.

With N = 30 and N = 120, and with cohorts of 1, the number of cohorts
changes (10, 40 and 60). `first_contour_cohort` counts cohorts, so it is not
comparable across those rows.

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
- `calibration_baseline.csv`, `calibration_alpha_0.10.csv`,
  `calibration_alpha_0.20.csv`, `calibration_shift_0.1.csv`,
  `calibration_shift_0.075.csv`, `calibration_shift_0.05.csv`,
  `calibration_shift_0.025.csv`, `riviere2014_scenario1_pcs.csv`.

## Seeds

Each configuration has its own `seed` in `mtd_configs()`. It runs its 10000
trials in one process, with `set.seed(seed)`. The parallel loop is over
configurations, so results do not depend on the number of cores or on which
other configurations run. The calibration uses seed 1.

## How to run

From this folder:

    Rscript calibrate.R          # baseline, alpha 0.05: about 10 minutes, one core
    Rscript calibrate.R 1 0.10 results/calibration_alpha_0.10.csv   # seed, alpha, file
    Rscript calibrate.R 1 0.05 results/calibration_shift_0.1.csv 0.1 0.1   # ..., band, shift
    Rscript run.R 4              # all configurations, about 4 minutes on 4 cores
    Rscript run.R 4 eps_0.5      # only these, and only their rows change
    Rscript report.R [id]        # tables; grids for one configuration

Run times for 10000 trials: about 25 seconds per row, `cohort_1` about 75
seconds, `N_120` about 50 seconds, `ref_ms` about 230 seconds (the safety
rule is slower). Each shift calibration takes about 10 minutes on one core; the four can run as
parallel processes. `run.R` needs the installed gPIPE.

`tests/testthat/test-simulation-mtd-sensitivity.R` checks the stored results.
