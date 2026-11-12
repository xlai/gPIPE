# Cheung (2022) Table 2: gPIPE results

Generated from the stored results in `results/`. Settings and the pass rule are in `README.md`. 400 trials per scenario (100,800 per design); Cheung used 40 per scenario.

## Summary

- **Pooled posterior (gPIPE's model, the pre-specified primary run): 6 of 15 cells pass.** All five FDR cells pass and so does n_eff for the non-adaptive design. Every TPR cell fails, by +0.15 for the non-adaptive design down to +0.03 for 2-PP, and the four adaptive n_eff cells fail by 1 to 4 patients. gPIPE is more powerful than the printed numbers at the same epsilon.
- **Per-cell posterior (Cheung's Sec. 2.2 model, a setting the paper leaves open; variant run after the primary): 6 of 15 cells pass.** The TPR gaps shrink to -0.03 to +0.01, but the FDR is now a little low (0.040 against 0.048 for 1-EQ and 1-PP) and most cells still miss the 3 SE margin, because the margin is tight (about 0.01 for TPR). The remaining differences are at most 0.03 in TPR, 0.009 in FDR and 2 patients in n_eff.
- The non-adaptive design (no allocation rule) already misses in TPR with the pooled posterior, so the gap comes from the posterior p_jk, not from the allocation rules or the FDR/TPR code. The script's FDR and TPR equal `summary()` on the same trials (test).
- Under both posteriors the TPR ranking of the five designs is Cheung's (non-adaptive < 1-PP < 1-EQ < 2-EQ < 2-PP).
- The margin was not widened. Following the plan in #31, the remaining gap is documented and agreement is claimed only qualitatively (design ranking, size of the effects).

## Pooled posterior (primary)

| design | metric | printed | gPIPE | diff | SE | z | 3 SE margin | pass |
|---|---|---|---|---|---|---|---|---|
| non-adaptive | FDR | 0.0460 | 0.0494 | 0.0034 | 0.0017 | 2.0 | 0.0050 | yes |
| non-adaptive | TPR | 0.680 | 0.830 | 0.150 | 0.004 | 42.5 | 0.011 | **no** |
| non-adaptive | n_eff | 63.0 | 62.5 | -0.5 | 0.3 | -1.7 | 0.9 | yes |
| 1-EQ | FDR | 0.0480 | 0.0451 | -0.0029 | 0.0016 | -1.9 | 0.0047 | yes |
| 1-EQ | TPR | 0.790 | 0.859 | 0.069 | 0.003 | 21.8 | 0.010 | **no** |
| 1-EQ | n_eff | 103.0 | 98.6 | -4.4 | 0.3 | -14.6 | 0.9 | **no** |
| 1-PP | FDR | 0.0480 | 0.0470 | -0.0010 | 0.0016 | -0.6 | 0.0048 | yes |
| 1-PP | TPR | 0.780 | 0.843 | 0.063 | 0.003 | 19.7 | 0.010 | **no** |
| 1-PP | n_eff | 111.0 | 108.3 | -2.7 | 0.3 | -9.0 | 0.9 | **no** |
| 2-EQ | FDR | 0.0460 | 0.0459 | -0.0001 | 0.0016 | -0.1 | 0.0047 | yes |
| 2-EQ | TPR | 0.880 | 0.933 | 0.053 | 0.003 | 17.1 | 0.009 | **no** |
| 2-EQ | n_eff | 66.0 | 67.0 | 1.0 | 0.3 | 3.3 | 0.9 | **no** |
| 2-PP | FDR | 0.0490 | 0.0514 | 0.0024 | 0.0017 | 1.4 | 0.0050 | yes |
| 2-PP | TPR | 0.910 | 0.945 | 0.035 | 0.003 | 11.4 | 0.009 | **no** |
| 2-PP | n_eff | 91.0 | 88.5 | -2.5 | 0.3 | -7.7 | 1.0 | **no** |

## Per-cell posterior (variant)

| design | metric | printed | gPIPE | diff | SE | z | 3 SE margin | pass |
|---|---|---|---|---|---|---|---|---|
| non-adaptive | FDR | 0.0460 | 0.0419 | -0.0041 | 0.0016 | -2.5 | 0.0049 | yes |
| non-adaptive | TPR | 0.680 | 0.650 | -0.030 | 0.004 | -7.3 | 0.012 | **no** |
| non-adaptive | n_eff | 63.0 | 62.5 | -0.5 | 0.3 | -1.7 | 0.9 | yes |
| 1-EQ | FDR | 0.0480 | 0.0396 | -0.0084 | 0.0015 | -5.5 | 0.0046 | **no** |
| 1-EQ | TPR | 0.790 | 0.769 | -0.021 | 0.003 | -6.6 | 0.010 | **no** |
| 1-EQ | n_eff | 103.0 | 104.4 | 1.4 | 0.3 | 4.6 | 0.9 | **no** |
| 1-PP | FDR | 0.0480 | 0.0395 | -0.0085 | 0.0015 | -5.5 | 0.0046 | **no** |
| 1-PP | TPR | 0.780 | 0.761 | -0.019 | 0.003 | -5.8 | 0.010 | **no** |
| 1-PP | n_eff | 111.0 | 111.8 | 0.8 | 0.3 | 2.8 | 0.9 | yes |
| 2-EQ | FDR | 0.0460 | 0.0437 | -0.0023 | 0.0016 | -1.5 | 0.0047 | yes |
| 2-EQ | TPR | 0.880 | 0.876 | -0.004 | 0.003 | -1.3 | 0.010 | yes |
| 2-EQ | n_eff | 66.0 | 68.1 | 2.1 | 0.3 | 6.6 | 0.9 | **no** |
| 2-PP | FDR | 0.0490 | 0.0488 | -0.0002 | 0.0016 | -0.1 | 0.0049 | yes |
| 2-PP | TPR | 0.910 | 0.920 | 0.010 | 0.003 | 3.1 | 0.009 | **no** |
| 2-PP | n_eff | 91.0 | 89.4 | -1.6 | 0.3 | -4.9 | 1.0 | **no** |

`SE` is the combined SE (Cheung's, gPIPE's and rounding). The 3 SE margin is the smallest difference the check can detect: about 0.005 for FDR, 0.010 for TPR and 0.9 patients for n_eff.

## Epsilon at FDR <= 0.05 (reported, not a gate)

Non-adaptive design, Cheung's grid 0.01 to 0.10 step 0.0025, 400 trials per scenario:

- pooled posterior: 0.0525 (printed 0.0525)
- per-cell posterior: 0.0575 (printed 0.0525)

Adaptive designs, `calibrate_eps()` (pooled posterior, bisection to 0.0025, 40 trials per scenario as Cheung; the FDR curve is noisy, so one grid step of 0.0025 or a bit more is expected):

| design | calibrated | printed | difference in grid steps | FDR | TPR |
|---|---|---|---|---|---|
| 1-EQ | 0.0536 | 0.0500 | 1.4 | 0.049 | 0.868 |
| 1-PP | 0.0648 | 0.0625 | 0.9 | 0.049 | 0.850 |
| 2-EQ | 0.0452 | 0.0425 | 1.1 | 0.050 | 0.937 |
| 2-PP | 0.0522 | 0.0525 | -0.1 | 0.048 | 0.943 |

