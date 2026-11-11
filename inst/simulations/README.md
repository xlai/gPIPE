# Simulation checks against published results

Scripts here rerun designs from the papers and compare gPIPE with the printed
numbers. Stored results are in `results/`; the tests in `tests/testthat` read
them and are skipped on CRAN.

## `cheung-table2.R`: Cheung, Chandereng and Diaz (2022), Table 2

Settings follow Sec. 4.1 of the paper: 5 x 5 grid, N = 125, cohort 1, 252
plateau scenarios (every up-set of the grid; mean 0.15 for an effective
combination, 0 for an ineffective one; true SD 0.16), `theta = 0`. Five
designs at the printed epsilon (non-adaptive, 1-EQ, 1-PP, 2-EQ, 2-PP). Cheung
used 40 trials per scenario; we use 400.

Pass rule (issue #31, step 7): a cell (design x FDR, TPR, n_eff) passes when
`|gPIPE - printed| <= 3 * sqrt(SE_Cheung^2 + SE_gPIPE^2 + SE_rounding^2)`.
All 15 cells must pass. The margin is not widened after seeing results.

Settings the paper leaves open, listed before the run (the only ones that may
be varied if a cell fails):

1. **Posterior for the normal outcome.** Cheung's sigma handling is not
   published. gPIPE uses a pooled-variance t posterior under a flat prior,
   with `sigma0 = 0.16` until a variance can be estimated. This is the primary
   run. The variant `cell` gives every combination its own variance (flat prior
   on the mean, 1/sigma on the SD, so a t posterior with n - 1 degrees of
   freedom; known SD below 2 patients), which is how Sec. 2.2 of the paper reads.
2. **TPR in the all-ineffective scenario** (eq. 8 is 0/0). Primary: that
   scenario is left out of the TPR average (as `summary()` does). Variants
   reported: count it as 0, or as 1.
3. **n_eff** is the mean number of patients treated at a truly effective
   combination, over all 252 scenarios (0 for the all-ineffective scenario).
4. **Non-adaptive design:** 5 patients on each of the 25 combinations, final
   estimate at the printed epsilon, weights n + 1 (constant, so equal weights).

Run:

    Rscript inst/simulations/cheung-table2.R main   [reps]   # five designs, pooled posterior
    Rscript inst/simulations/cheung-table2.R main 400 cell   # same, per-cell posterior (setting 1)
    Rscript inst/simulations/cheung-table2.R calib           # soft epsilon check

The calibrated epsilon is reported only; it is not a pass/fail gate. Results, with
the two 15-cell tables, are in `cheung-table2-results.md`.
