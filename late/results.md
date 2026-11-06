# Nested ANOVA: run > well > electrode, late DIV [16, 17, 18], vehicle control, y = ln(mean ISI [ms])

Thinning: min pairwise distance d = 150 um; wells kept with m >= 10 thinned electrodes: 42 of 42 wells, 7 runs, N = 4454 electrodes. Unbalanced design (m per well 38-202).

## ANOVA table (thinned); p = parametric F-test

| source            |        SS |   df |     MS |        F |        p | den          |
|:------------------|----------:|-----:|-------:|---------:|---------:|:-------------|
| Run               |    7.0545 |    6 | 1.1757 |   1.7198 |   0.1454 | MS_well(run) |
| Well(Run)         |   23.9281 |   35 | 0.6837 |   2.3775 |   0      | MS_electrode |
| Electrode (resid) | 1268.67   | 4412 | 0.2875 | nan      | nan      |              |

Permutation p (5000 perms): run = 0.1534, well(run) = 0.0004

## Per-run well test (thinned; parametric F with pooled MS_electrode)

| run   |   wells |      F |   df1 |      p |
|:------|--------:|-------:|------:|-------:|
| 22    |       6 | 1.6999 |     5 | 0.131  |
| 23    |       6 | 3.1245 |     5 | 0.0081 |
| 24    |       6 | 0.7724 |     5 | 0.5695 |
| 42    |       6 | 2.1966 |     5 | 0.0519 |
| 43    |       6 | 2.5409 |     5 | 0.0265 |
| 45    |       6 | 5.2208 |     5 | 0.0001 |
| AR    |       6 | 1.0877 |     5 | 0.3649 |

## Variance components

| method               |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |     N |
|:---------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|------:|
| ANOVA thinned        |    0.0007 |     0.0038 |          0.2875 |     0.0025 |      0.013  |           0.9845 |  4454 |
| ANOVA all electrodes |    0.0026 |     0.0049 |          0.2123 |     0.0119 |      0.0224 |           0.9657 | 52658 |
| REML thinned         |    0.0007 |     0.004  |          0.2876 |     0.0023 |      0.0136 |           0.9841 |  4454 |
| REML all electrodes  |    0.0009 |     0.0076 |          0.2123 |     0.0042 |      0.0345 |           0.9613 | 52658 |