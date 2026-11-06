# Nested Frechet (W2) ANOVA of per-electrode ISI distributions — mid DIV [9, 10, 11], vehicle control

Per-electrode distribution: 2-component GMM on log10(ISI ms), converted to quantile function; electrodes need >= 30 ISIs; thinned to >= 150 um; wells kept with >= 10 thinned electrodes: 32 wells in 7 runs, N = 1281 electrodes (m per well 10-104). Quantile grid K = 200.

## Frechet ANOVA table (SS in (log10 ISI)^2 units, L2 on quantile functions)

| source            |       SS |   df |      MS |    F_stat |   perm_p |
|:------------------|---------:|-----:|--------:|----------:|---------:|
| Run               |  16.0533 |    6 | 2.67555 |   4.38899 |   0.0008 |
| Well(Run)         |  15.2401 |   25 | 0.6096  |   2.92138 |   0.0002 |
| Electrode (resid) | 260.628  | 1249 | 0.20867 | nan       | nan      |

F_stat is the usual ratio MS_run/MS_well(run) and MS_well(run)/MS_e; its null distribution is NOT F, p-values are from 5000 permutations.

## Per-run well test (permutation; eta2 = SS_well / (SS_well + SS_e) within run, Holm across runs)

| run   |   wells |   SS_well |      eta2 |         p |    p_holm |
|:------|--------:|----------:|----------:|----------:|----------:|
| 22    |       6 |   3.23269 |   0.17775 |   0.0002  |   0.001   |
| 23    |       1 | nan       | nan       | nan       | nan       |
| 24    |       1 | nan       | nan       | nan       | nan       |
| 42    |       6 |   3.07202 |   0.08852 |   0.002   |   0.006   |
| 43    |       6 |   1.29674 |   0.02109 |   0.29634 |   0.29634 |
| 45    |       6 |   5.76401 |   0.07364 |   0.0002  |   0.001   |
| AR    |       6 |   1.87464 |   0.0231  |   0.09158 |   0.18316 |

## Variance components (Frechet variance, (log10 ISI)^2) and share of total

| method                                             |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |    N |
|:---------------------------------------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|-----:|
| thinned (tested)                                   |   0.01258 |    0.00993 |         0.20867 |    0.05443 |     0.04295 |          0.90262 | 1281 |
| all eligible electrodes, empirical Q (descriptive) |   0.00761 |    0.01543 |         0.19081 |    0.03558 |     0.07214 |          0.89228 | 7888 |