# Nested Frechet (W2) ANOVA of per-electrode ISI distributions — early DIV [3, 4], vehicle control

Per-electrode distribution: 2-component GMM on log10(ISI ms), converted to quantile function; electrodes need >= 30 ISIs; thinned to >= 100 um; wells kept with >= 10 thinned electrodes: 16 wells in 6 runs, N = 364 electrodes (m per well 11-65). Quantile grid K = 200.

## Frechet ANOVA table (SS in (log10 ISI)^2 units, L2 on quantile functions)

| source            |       SS |   df |      MS |    F_stat |    perm_p |
|:------------------|---------:|-----:|--------:|----------:|----------:|
| Run               |  9.58583 |    5 | 1.91717 |   1.94715 |   0.09338 |
| Well(Run)         |  9.84601 |   10 | 0.9846  |   6.22505 |   0.0002  |
| Electrode (resid) | 55.0423  |  348 | 0.15817 | nan       | nan       |

F_stat is the usual ratio MS_run/MS_well(run) and MS_well(run)/MS_e; its null distribution is NOT F, p-values are from 5000 permutations.

## Per-run well test (permutation; eta2 = SS_well / (SS_well + SS_e) within run, Holm across runs)

| run   |   wells |   SS_well |      eta2 |         p |    p_holm |
|:------|--------:|----------:|----------:|----------:|----------:|
| 22    |       3 |   3.18383 |   0.27882 |   0.0002  |   0.001   |
| 23    |       1 | nan       | nan       | nan       | nan       |
| 42    |       2 |   1.63355 |   0.22873 |   0.0084  |   0.0168  |
| 43    |       2 |   2.03543 |   0.35118 |   0.0002  |   0.001   |
| 45    |       2 |   0.25584 |   0.03749 |   0.19876 |   0.19876 |
| AR    |       6 |   2.73736 |   0.08176 |   0.0036  |   0.0108  |

## Variance components (Frechet variance, (log10 ISI)^2) and share of total

| method                                             |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |   N |
|:---------------------------------------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|----:|
| thinned (tested)                                   |   0.02077 |    0.03749 |         0.15817 |    0.09599 |     0.17321 |          0.7308  | 364 |
| all eligible electrodes, empirical Q (descriptive) |   0.02382 |    0.02803 |         0.16791 |    0.1084  |     0.12754 |          0.76405 | 847 |