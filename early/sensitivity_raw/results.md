# Nested ANOVA: run > well > electrode, early DIV [3, 4], vehicle control, y = mean ISI [ms] (raw)

Thinning: min pairwise distance d = 100 um; wells kept with m >= 10 thinned electrodes: 37 of 41 wells, 7 runs, N = 1541 electrodes. Unbalanced design (m per well 10-116).

## ANOVA table (thinned); p = parametric F-test

| source            |              SS |   df |      MS |        F |        p | den          |
|:------------------|----------------:|-----:|--------:|---------:|---------:|:-------------|
| Run               | 21629.7         |    6 | 3604.94 |   1.3976 |   0.248  | MS_well(run) |
| Well(Run)         | 77380.5         |   30 | 2579.35 |   2.0371 |   0.0008 | MS_electrode |
| Electrode (resid) |     1.90435e+06 | 1504 | 1266.19 | nan      | nan      |              |

Permutation p (5000 perms): run = 0.2601, well(run) = 0.0030

## Per-run well test (thinned; parametric F with pooled MS_electrode)

| run   |   wells |      F |   df1 |      p |
|:------|--------:|-------:|------:|-------:|
| 22    |       6 | 2.1452 |     5 | 0.0577 |
| 23    |       6 | 3.8767 |     5 | 0.0017 |
| 24    |       3 | 3.3055 |     2 | 0.0369 |
| 42    |       6 | 1.9077 |     5 | 0.0901 |
| 43    |       5 | 2.0935 |     4 | 0.0794 |
| 45    |       5 | 1.2493 |     4 | 0.2881 |
| AR    |       6 | 0.2965 |     5 | 0.915  |

## Variance components

| method               |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |    N |
|:---------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|-----:|
| ANOVA thinned        |    4.5454 |    32.1308 |         1266.19 |     0.0035 |      0.0247 |           0.9718 | 1541 |
| ANOVA all electrodes |   11.8921 |    21.1116 |         1004.24 |     0.0115 |      0.0204 |           0.9682 | 5981 |
| REML thinned         |   36.1599 |   144.607  |         1248.95 |     0.0253 |      0.1011 |           0.8736 | 1541 |
| REML all electrodes  |  621.86   |    41.607  |         1004.03 |     0.3729 |      0.025  |           0.6021 | 5981 |