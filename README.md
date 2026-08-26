# Few-Treated Survival Simulation

This project implements simulation code for few-treated matching inference.

The current code covers three parts:

1. Ferman 2021 non-survival replication with nearest-neighbor matching and Ferman-style sign-change p-values.
2. Baba-Yoshida 2024 survival baseline with CEM matching and Gaussian weighted log-rank p-values.
3. Proposed survival Version A and Version B, using treated-level survival contributions and Ferman-style sign-change p-values.

## Main R Files

- `R/few_treated_survival_methods.R`: main implementation file. It is kept as one file, but organized into clear `MODULE` blocks for p-value construction, matching, survival methods, toy examples, simulation grids, and Ferman replication.
- `R/run_simulation.R`: simulation tasks, grid settings, and result-writing code.
- `R/test_toy_examples.R`: unit tests for toy examples and key function wrappers.
- `R/test_nn_match_diagnostics.R`: quick diagnostic test for nearest-neighbor reuse calculations.
- `R/test_cem_package_validation.R`: quick comparison between manual CEM grouping and `MatchIt::matchit(method = "cem")`.

## Main Commands

Run tests:

```r
system("Rscript R/run_simulation.R tests")
```

Run only the toy/unit tests:

```r
system("Rscript R/run_simulation.R unit_tests")
```

Print the full run plan:

```r
system("Rscript R/run_simulation.R plan")
```

Run the full Ferman, Baba-Yoshida, and Version A/B simulation. The Ferman part now includes both null and alternative settings:

```r
system("Rscript R/run_simulation.R complete_300")
```

Run only the Ferman null-plus-alternative grid:

```r
system("Rscript R/run_simulation.R ferman_300")
```

Run the Ferman Table-2-aligned check:

```r
system("Rscript R/run_simulation.R ferman_table2_300")
```

This uses `alpha = 0.10` because Ferman 2021 Table 2 reports 10% rejection rates. The regular `ferman_300` task is kept at `alpha = 0.05` as the main project run. The Table 2 comparison should use the deduplicated null rows written by the `ferman_table2_300` task.

Run only the old Ferman null-only grid:

```r
system("Rscript R/run_simulation.R ferman_null_300")
```

Run the Ferman covariate-dimension extension. The pilot is a faster check on Panel A; the full version runs Panels A-E with dimensions 1, 2, 4, and 6:

```r
system("Rscript R/run_simulation.R ferman_dim_pilot 50")
system("Rscript R/run_simulation.R ferman_dim_300")
```

Run only the Version A/B 30-minute pilot:

```r
system("Rscript R/run_simulation.R version_ab_30min_pilot")
```

Run only Version A with the proposed survival grid:

```r
system("Rscript R/run_simulation.R version_a_300_parallel")
```

## Current Notes

- CEM matching is implemented through `MatchIt::matchit(method = "cem")`.
- CEM bins use `max(4, floor((N0 + N1)^(1/3)))`.
- Ferman alternative is controlled by an additive treated-outcome effect `tau = 0.5`; null uses `tau = 0`.
- The main Ferman run uses `alpha = 0.05`. For direct comparison with Ferman 2021 Table 2, use `ferman_table2_300`, which uses `alpha = 0.10` and writes a null-only Table 2 comparison file.
- The Ferman covariate-dimension extension keeps the original Panel A-E outcome model tied to the first covariate and adds extra matching covariates. This checks whether nearest-neighbor matching degrades as matching dimension increases, without changing the original one-dimensional outcome surface.
- Other Ferman tables are not forced to align unless requested separately; Table 2 is the main comparison target for the current nearest-neighbor few-treated replication.
- When `N1` is very small, sign-change p-values are discrete. With four or fewer effective sign units, the 5% critical value is the maximum randomization statistic, so rejection can be exactly zero.
- Baba-Yoshida baseline uses Gaussian p-values, so sign-change diagnostics such as `mean_sign_units` and `mean_contributions` are not applicable there.
- Version A/B full-grid simulation is computationally expensive on a laptop. The current package includes a partial full-grid result and a smaller pilot result for trend checking.

## Results To Review

Current selected result files are under:

```text
results/latest_results/
```

Grid-level checkpoint/latest files are under:

```text
results/grid_level/
```

Every new CSV output is also copied into a timestamped run folder:

```text
results/runs/<YYYYMMDD_HHMMSS>/
```

Old scattered result files are archived under:

```text
results/archive/
```
