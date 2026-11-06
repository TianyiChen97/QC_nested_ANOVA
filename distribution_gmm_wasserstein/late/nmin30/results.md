# Nested Frechet (W2) ANOVA of per-electrode ISI distributions — late DIV [16, 17, 18], vehicle control

Per-electrode distribution: 2-component GMM on log10(ISI ms), converted to quantile function; electrodes need >= 30 ISIs; thinned to >= 150 um; wells kept with >= 10 thinned electrodes: 41 wells in 7 runs, N = 2181 electrodes (m per well 13-144). Quantile grid K = 200.

## Frechet ANOVA table (SS in (log10 ISI)^2 units, L2 on quantile functions)

| source            |       SS |   df |      MS |    F_stat |   perm_p |
|:------------------|---------:|-----:|--------:|----------:|---------:|
| Run               |  18.2364 |    6 | 3.0394  |   3.34865 |   0.0248 |
| Well(Run)         |  30.8601 |   34 | 0.90765 |   4.70989 |   0.0002 |
| Electrode (resid) | 412.402  | 2140 | 0.19271 | nan       | nan      |

F_stat is the usual ratio MS_run/MS_well(run) and MS_well(run)/MS_e; its null distribution is NOT F, p-values are from 5000 permutations.

## Per-run well test (permutation; eta2 = SS_well / (SS_well + SS_e) within run, Holm across runs)

| run   |   wells |   SS_well |    eta2 |       p |   p_holm |
|:------|--------:|----------:|--------:|--------:|---------:|
| 22    |       6 |   1.00958 | 0.04006 | 0.15337 |  0.21596 |
| 23    |       6 |   4.00811 | 0.13018 | 0.0002  |  0.0014  |
| 24    |       5 |   1.56373 | 0.08453 | 0.0246  |  0.07379 |
| 42    |       6 |   1.68497 | 0.04728 | 0.10798 |  0.21596 |
| 43    |       6 |   6.89058 | 0.05849 | 0.0002  |  0.0014  |
| 45    |       6 |  12.6951  | 0.13734 | 0.0002  |  0.0014  |
| AR    |       6 |   3.00801 | 0.02447 | 0.0024  |  0.0096  |

## Variance components (Frechet variance, (log10 ISI)^2) and share of total

| method                                             |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |     N |
|:---------------------------------------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|------:|
| thinned (tested)                                   |   0.00695 |    0.01388 |         0.19271 |    0.03254 |     0.06499 |          0.90248 |  2181 |
| all eligible electrodes, empirical Q (descriptive) |   0.00784 |    0.00806 |         0.15829 |    0.04501 |     0.04624 |          0.90874 | 19134 |