# Mander & Sweeting (2015), Tables III and IV

Validation of gPIPE's original PIPE design (`ms_design()`) against the
authors' package, pipe.design 0.5.1, and against the printed tables
(Statistics in Medicine 34:1261, <doi:10.1002/sim.6434>).

- `common.R`: scenarios (Tables I and II), the 16 + 7 configurations, and the
  summary both runs share.
- `run_gpipe.R`, `run_pipe_design.R`: the two runs, with 4000 trials per
  Table III configuration and 8000 per Table IV scenario (4 times the paper).
- `compare.R`: differences in Monte Carlo standard errors.
- `results/`: the stored runs and the printed tables (`ms2015_paper.csv`).

Design: weak prior (strength 1 / K per dose; medians from scenario 1 or A),
neighbouring constraint, safety rule at 0.8, start at the lowest
combination. Table III: 6 x 6, theta = 0.3, cohorts of 2, N = 40; closest
or adjacent doses; smallest sample size (`min_S`) or randomisation weighted
by 1 / S (`inverse_S`). Table IV: 4 x 4, theta = 0.2, cohorts of 1, N = 50,
closest, `min_S`.

Metrics, as in the paper: experimentation and recommendation percentages
by true-toxicity band, pooled over trials, and the mean number of
recommended doses. In Table IV a stopped trial counts once as "None" in the
recommendations, and its untreated patients count as "None" in the
experimentation.

Run from this folder (gPIPE about 25 minutes, pipe.design about 40 minutes
on 4 cores):

    Rscript run_gpipe.R 4 4
    Rscript run_pipe_design.R 4 4
    Rscript compare.R

`tests/testthat/test-validation-ms2015.R` checks the stored results.
