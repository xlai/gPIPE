# Golden-test fixtures (#17)

Known answers that do not come from gPIPE. The tests that use them run once
per implementation (see `../helper-impl.R`).

## `pipe_design_V5.json`, `pipe_design_V6.json`

Single decision states from pipe.design 0.5.1 (the Mander & Sweeting 2015
reference implementation), recorded with its own internals (`mtc.create`,
`monotonic.matrices`, `beta.med`).

- V5: 6 x 6 grid, theta = 0.3, cohort 2, N = 40, weak prior (a + b = 1/36),
  scenario-1 medians; 30 states. Reproduces Table III.
- V6: 4 x 4 grid, theta = 0.2, cohort 1, N = 50, weak prior (a + b = 1/16),
  scenario-A medians; 33 states. Reproduces Table IV.
- Safety threshold 0.8, neighbouring constraint, closest or adjacent rule.
- States come from real pipe.design runs plus six hand-made edge cases per
  setting (`edge-*`).

Indexing is 0-based `[doseA][doseB]`. `p` is P(pi <= theta | data), the
complement of gPIPE's p. Use the prior `a`, `b` stored in `meta`, not a
re-solved prior: at doses whose prior median equals theta, pipe.design's
solver leaves p about 1e-4 away from 0.5, which decides the mode contour.

pipe.design behaviours the tests pin down:

- The trial stops when no dose in the neighbourhood is safe, even if a safe
  dose exists elsewhere. pipe.design's nearest-safe-dose fallback never fires.
- Min-S ties are ties in n: pipe.design compares 1/S exactly, while a + b
  differs between doses by rounding only.

Generator, spec and the full Monte Carlo comparison with the paper live in the
project's validation folder (`pipe_design_reference/`, `make_snapshots.R`).

## `cheung2022_table1.csv`

Cheung, Chandereng & Diaz (2022), Table 1: 12 subjects of the sedentary-break
example on a 5 x 5 grid, for Algorithm 2-EQ (taper to eps = 0.0425, N = 125)
and Algorithm 1-EQ (eps = 0.05). Levels are 1-based.

- `p_printed`: the posterior probability printed for the dose just treated.
- `p_used`: the same, except 1-EQ step 2. The paper prints 0.95, which
  equals 1 - eps, so the estimate there sits on a knife edge. The value
  0.9505 is our substitution for the rounded print, not the paper's; it
  reproduces the paper's MED set.
- `med_hat`: the printed MED set after that subject, `;`-separated.

With w = n + 1 and untried p = 0.5, all 24 MED sets reproduce, and every dose
lies in S1 u S0 (Algorithm 1, step 3) of the previous estimate.
