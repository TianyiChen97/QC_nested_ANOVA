# Nested ANOVA: run > well > electrode, late DIV [16, 17, 18], vehicle control, y = mean ISI [ms] (raw)

Thinning: min pairwise distance d = 150 um; wells kept with m >= 10 thinned electrodes: 42 of 42 wells, 7 runs, N = 4454 electrodes. Unbalanced design (m per well 38-202).

## ANOVA table (thinned); p = parametric F-test

| source            |               SS |   df |      MS |        F |        p | den          |
|:------------------|-----------------:|-----:|--------:|---------:|---------:|:-------------|
| Run               |  30607           |    6 | 5101.16 |   1.2284 |   0.3155 | MS_well(run) |
| Well(Run)         | 145346           |   35 | 4152.75 |   3.1051 |   0      | MS_electrode |
| Electrode (resid) |      5.90055e+06 | 4412 | 1337.39 | nan      | nan      |              |

Permutation p (5000 perms): run = 0.3187, well(run) = 0.0002

## Per-run well test (thinned; parametric F with pooled MS_electrode)

| run   |   wells |      F |   df1 |      p |
|:------|--------:|-------:|------:|-------:|
| 22    |       6 | 2.3866 |     5 | 0.0359 |
| 23    |       6 | 3.8616 |     5 | 0.0017 |
| 24    |       6 | 1.3193 |     5 | 0.2526 |
| 42    |       6 | 2.3813 |     5 | 0.0363 |
| 43    |       6 | 3.8039 |     5 | 0.0019 |
| 45    |       6 | 6.6369 |     5 | 0      |
| AR    |       6 | 1.3462 |     5 | 0.2417 |

## Variance components

| method               |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |     N |
|:---------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|------:|
| ANOVA thinned        |    1.1295 |    26.9697 |         1337.39 |     0.0008 |      0.0198 |           0.9794 |  4454 |
| ANOVA all electrodes |   13.7975 |    25.2904 |         1023.71 |     0.013  |      0.0238 |           0.9632 | 52658 |
| REML thinned         |    1.1396 |    28.4785 |         1337.71 |     0.0008 |      0.0208 |           0.9783 |  4454 |
| REML all electrodes  |   11.0936 |    37.0984 |         1023.72 |     0.0103 |      0.0346 |           0.955  | 52658 |