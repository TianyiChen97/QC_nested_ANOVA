# Nested ANOVA: run > well > electrode, early DIV [3, 4], vehicle control, y = ln(mean ISI [ms])

Thinning: min pairwise distance d = 100 um; wells kept with m >= 10 thinned electrodes: 37 of 41 wells, 7 runs, N = 1541 electrodes. Unbalanced design (m per well 10-116).

## ANOVA table (thinned); p = parametric F-test

| source            |       SS |   df |     MS |        F |        p | den          |
|:------------------|---------:|-----:|-------:|---------:|---------:|:-------------|
| Run               |   9.8462 |    6 | 1.641  |   1.0479 |   0.4149 | MS_well(run) |
| Well(Run)         |  46.9822 |   30 | 1.5661 |   3.747  |   0      | MS_electrode |
| Electrode (resid) | 628.6    | 1504 | 0.418  | nan      | nan      |              |

Permutation p (5000 perms): run = 0.2803, well(run) = 0.0002

## Per-run well test (thinned; parametric F with pooled MS_electrode)

| run   |   wells |       F |   df1 |      p |
|:------|--------:|--------:|------:|-------:|
| 22    |       6 |  4.5719 |     5 | 0.0004 |
| 23    |       6 | 10.9326 |     5 | 0      |
| 24    |       3 |  8.1423 |     2 | 0.0003 |
| 42    |       6 |  0.7758 |     5 | 0.5671 |
| 43    |       5 |  1.1241 |     4 | 0.3434 |
| 45    |       5 |  0.5326 |     4 | 0.7118 |
| AR    |       6 |  1.6195 |     5 | 0.1517 |

## Variance components

| method               |   var_run |   var_well |   var_electrode |   prop_run |   prop_well |   prop_electrode |    N |
|:---------------------|----------:|-----------:|----------------:|-----------:|------------:|-----------------:|-----:|
| ANOVA thinned        |    0      |     0.0281 |          0.418  |     0.0001 |      0.063  |           0.9369 | 1541 |
| ANOVA all electrodes |    0.0012 |     0.0079 |          0.2295 |     0.005  |      0.0331 |           0.9619 | 5981 |
| REML thinned         |    0.0053 |     0.0318 |          0.4191 |     0.0117 |      0.0697 |           0.9187 | 1541 |
| REML all electrodes  |    0.0019 |     0.0205 |          0.2297 |     0.0075 |      0.0814 |           0.9111 | 5981 |