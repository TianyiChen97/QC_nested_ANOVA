# Nested Frechet (2-Wasserstein) ANOVA of per-electrode ISI distributions

Each electrode is a distribution, not one number: a 2-component Gaussian mixture is fitted to log10(ISI ms) (5 parameters: w, mu1, sigma1, mu2, sigma2), converted to its quantile function Q(u) on a grid of 200 levels.
In 1D, W2^2(P1,P2) = ||Q1 - Q2||^2_L2, so the data live in a Hilbert space: the barycenter of a set of distributions = pointwise mean of their quantile functions
(invariant to mixture label switching), and SS_total = SS_run + SS_well(run) + SS_electrode holds exactly, with the same Searle coefficients as the scalar analysis.
Tests are permutation tests only (no F reference for functional data): run = permute run labels of whole wells; well = permute electrodes among wells within run; per-run well tests with Holm.

Folders `early/ mid/ late/`, each with `nmin30/` (electrodes need >= 30 ISIs) containing:
`frechet_tree_<tp>_vc.png` (main figure), `results.md`, `frechet_anova_table.csv`, `per_run_well_tests.csv`, `variance_components.csv`, `variance_profile_by_quantile_level.csv`,
`electrode_fits.csv` (GMM parameters per electrode + BIC + W2 fit-vs-data), `Q_gmm.npy`/`Q_emp.npy`, `well_counts.csv`, and `sensitivity_empirical_quantiles/` (same analysis with each electrode's EMPIRICAL quantile function, no mixture model).

Code (`code/`, run in order with early|mid|late): 01_extract_isis.py -> 02_fit_gmm_quantiles.py <tp> [N_MIN] -> 03_frechet_nested_anova.py <tp> [N_MIN] [gmm|emp] -> 04_tree_figure.py <tp> [N_MIN] [gmm|emp].
Thinning distances follow the residual correlograms of the mean-ISI analysis (early 100 um, mid/late 150 um); thinning uses positions only among eligible (>= 30 ISIs) electrodes.

## Post hoc (code/06_posthoc.py, output `<tp>/nmin30/posthoc/`)
Max-T permutation tests (5000 perms) that control the family-wise error within each family.
Wells: pairwise and one-vs-rest within each run (family = wells of that run; electrodes permuted among the run's wells). Runs: pairwise and one-vs-rest (family = all runs; units = wells; run labels of whole wells permuted).
Files: well_pairs.csv, well_vs_rest_of_run.csv, run_pairs.csv, run_vs_rest.csv, posthoc_<tp>.png. W2_distance is in log10 ISI units.
Run-level contrasts use full-range ISIs and inherit the 30 s vs 60 s scan-length confound; well-level contrasts do not.
