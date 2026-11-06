# Nested ANOVA: run > well > electrode, mid DIV [9, 10, 11], vehicle control, y = mean ISI [ms] (raw)

Thinning: min pairwise distance d = 150 um; wells kept with m >= 10 thinned electrodes: 41 of 42 wells, 7 runs, N = 2919 electrodes. Unbalanced design (m per well 26-147).

## ANOVA table (thinned); p = parametric F-test

| source            |              SS |   df |      MS |        F |        p | den          |
|:------------------|----------------:|-----:|--------:|---------:|---------:|:-------------|
| Run               | 56756.3         |    6 | 9459.39 |   5.0445 |   0.0009 | MS_well(run) |
| Well(Run)         | 63756           |   34 | 1875.18 |   1.4056 |   0.0601 | MS_electrode |
| Electrode (resid) |     3.83957e+06 | 2878 | 1334.11 | nan      | nan      |              |

Permutation p (5000 perms): run = 0.0024, well(run) = 0.0744

## Per-run well test (thinned; parametric F with pooled MS_electrode)

| run   |   wells |      F |   df1 |      p |
|:------|--------:|-------:|------:|-------:|
| 22    |       6 | 2.1082 |     5 | 0.0616 |
| 23    |       6 | 3.4089 |     5 | 0.0045 |
| 24    |       5 | 0.2903 |     4 | 0.8844 |
| 42    |       6 | 1.5357 |     5 | 0.1752 |
| 43    |       6 | 0.3857 |     5 | 0.8589 |
| 45    |       6 | 0.8027 |     5 | 0.5476 |
| AR    |       6 | 1.0843 |     5 | 0.3669 |

## Variance components

| method               |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |     N |
|:---------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|------:|
| ANOVA thinned        |   18.3586 |     7.6975 |         1334.11 |     0.0135 |      0.0057 |           0.9808 |  2919 |
| ANOVA all electrodes |   34.8493 |    10.223  |         1118.69 |     0.0299 |      0.0088 |           0.9613 | 25989 |
| REML thinned         |   15.433  |     6.8031 |         1334.93 |     0.0114 |      0.005  |           0.9836 |  2919 |
| REML all electrodes  |   16.1056 |    17.4016 |         1118.83 |     0.014  |      0.0151 |           0.9709 | 25989 |