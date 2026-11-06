# Nested Frechet (W2) ANOVA of per-electrode ISI distributions — late DIV [16, 17, 18], vehicle control

Per-electrode distribution: EMPIRICAL quantile function of log10(ISI ms) (sensitivity: no mixture model); electrodes need >= 30 ISIs; thinned to >= 150 um; wells kept with >= 10 thinned electrodes: 41 wells in 7 runs, N = 2181 electrodes (m per well 13-144). Quantile grid K = 200.

## Frechet ANOVA table (SS in (log10 ISI)^2 units, L2 on quantile functions)

| source            |       SS |   df |      MS |    F_stat |   perm_p |
|:------------------|---------:|-----:|--------:|----------:|---------:|
| Run               |  19.2442 |    6 | 3.20736 |   3.4828  |   0.019  |
| Well(Run)         |  31.3111 |   34 | 0.92091 |   4.70513 |   0.0002 |
| Electrode (resid) | 418.853  | 2140 | 0.19573 | nan       | nan      |

F_stat is the usual ratio MS_run/MS_well(run) and MS_well(run)/MS_e; its null distribution is NOT F, p-values are from 5000 permutations.

## Per-run well test (permutation; eta2 = SS_well / (SS_well + SS_e) within run, Holm across runs)

| run   |   wells |   SS_well |    eta2 |       p |   p_holm |
|:------|--------:|----------:|--------:|--------:|---------:|
| 22    |       6 |   0.91537 | 0.03654 | 0.21276 |  0.21596 |
| 23    |       6 |   4.15518 | 0.13133 | 0.0002  |  0.0014  |
| 24    |       5 |   1.61115 | 0.08601 | 0.0244  |  0.07319 |
| 42    |       6 |   1.71603 | 0.04748 | 0.10798 |  0.21596 |
| 43    |       6 |   6.95168 | 0.05784 | 0.0002  |  0.0014  |
| 45    |       6 |  12.9816  | 0.13959 | 0.0002  |  0.0014  |
| AR    |       6 |   2.98009 | 0.02376 | 0.0036  |  0.0144  |

## Variance components (Frechet variance, (log10 ISI)^2) and share of total

| method                                             |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |     N |
|:---------------------------------------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|------:|
| thinned (tested)                                   |   0.00747 |    0.01408 |         0.19573 |    0.03439 |     0.06478 |          0.90083 |  2181 |
| all eligible electrodes, empirical Q (descriptive) |   0.00784 |    0.00806 |         0.15829 |    0.04501 |     0.04624 |          0.90874 | 19134 |