# Few-Treated Survival Simulation

This project implements simulation code for few-treated matching inference.

The current code covers three parts:

1. Ferman 2021 non-survival replication with nearest-neighbor matching and Ferman-style sign-change p-values.
2. Baba-Yoshida 2024 survival baseline with CEM matching and Gaussian weighted log-rank p-values.
3. Proposed survival Version A and Version B, using treated-level survival contributions and Ferman-style sign-change p-values.

## Main R Files

- `R/few_treated_survival_methods.R`: method implementations, matching, p-value construction, diagnostics, and data-generating processes.
- `R/run_simulation.R`: simulation tasks, grid settings, and result-writing code.
- `R/test_nn_match_diagnostics.R`: quick diagnostic test for nearest-neighbor reuse calculations.
- `R/test_cem_package_validation.R`: quick comparison between manual CEM grouping and `MatchIt::matchit(method = "cem")`.

## Main Commands

Run tests:

```r
system("Rscript R/run_simulation.R tests")
```

Print the full run plan:

```r
system("Rscript R/run_simulation.R plan")
```

Run the full Ferman, Baba-Yoshida, and Version A/B simulation:

```r
system("Rscript R/run_simulation.R complete_300")
```

Run only the Version A/B 30-minute pilot:

```r
system("Rscript R/run_simulation.R version_ab_30min_pilot")
```

## Current Notes

- CEM matching is implemented through `MatchIt::matchit(method = "cem")`.
- CEM bins use `max(4, floor((N0 + N1)^(1/3)))`.
- Baba-Yoshida baseline uses Gaussian p-values, so sign-change diagnostics such as `mean_sign_units` and `mean_contributions` are not applicable there.
- Version A/B full-grid simulation is computationally expensive on a laptop. The current package includes a partial full-grid result and a smaller pilot result for trend checking.

## Results To Review

Selected result files are under:

```text
results/review_package/00_files_to_send_teacher/
```
