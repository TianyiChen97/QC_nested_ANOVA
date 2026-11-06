# Nested ANOVA: run > well > electrode, mid DIV [9, 10, 11], vehicle control, y = ln(mean ISI [ms])

Thinning: min pairwise distance d = 150 um; wells kept with m >= 10 thinned electrodes: 41 of 42 wells, 7 runs, N = 2919 electrodes. Unbalanced design (m per well 26-147).

## ANOVA table (thinned); p = parametric F-test

| source            |       SS |   df |     MS |        F |        p | den          |
|:------------------|---------:|-----:|-------:|---------:|---------:|:-------------|
| Run               |  14.9282 |    6 | 2.488  |   3.5147 |   0.0082 | MS_well(run) |
| Well(Run)         |  24.0681 |   34 | 0.7079 |   2.287  |   0      | MS_electrode |
| Electrode (resid) | 890.81   | 2878 | 0.3095 | nan      | nan      |              |

Permutation p (5000 perms): run = 0.0062, well(run) = 0.0006

## Per-run well test (thinned; parametric F with pooled MS_electrode)

| run   |   wells |      F |   df1 |      p |
|:------|--------:|-------:|------:|-------:|
| 22    |       6 | 1.7357 |     5 | 0.123  |
| 23    |       6 | 9.1891 |     5 | 0      |
| 24    |       5 | 0.2408 |     4 | 0.9153 |
| 42    |       6 | 0.9355 |     5 | 0.4567 |
| 43    |       6 | 0.4457 |     5 | 0.8166 |
| 45    |       6 | 0.8538 |     5 | 0.5115 |
| AR    |       6 | 2.1993 |     5 | 0.0518 |

## Variance components

| method               |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |     N |
|:---------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|------:|
| ANOVA thinned        |    0.0043 |     0.0057 |          0.3095 |     0.0134 |      0.0177 |           0.9689 |  2919 |
| ANOVA all electrodes |    0.0066 |     0.0017 |          0.2297 |     0.0278 |      0.007  |           0.9651 | 25989 |
| REML thinned         |    0.004  |     0.0057 |          0.3099 |     0.0125 |      0.0179 |           0.9696 |  2919 |
| REML all electrodes  |    0.0037 |     0.0032 |          0.2297 |     0.0158 |      0.0135 |           0.9707 | 25989 |