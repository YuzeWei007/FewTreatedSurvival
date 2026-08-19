## Quick check for nn_match_diagnostics()
##
## Run from either the project folder or the R/ folder:
## Rscript R/test_nn_match_diagnostics.R

if (file.exists(file.path("R", "few_treated_survival_methods.R"))) {
  source(file.path("R", "few_treated_survival_methods.R"))
} else {
  source("few_treated_survival_methods.R")
}

dat <- toy_scalar_data()

cat("Toy data:\n")
print(dat)

matches_m1 <- nn_match(
  data = dat,
  covariates = "X",
  M = 1,
  scaling = FALSE
)

cat("\nNearest-neighbor matches, M = 1:\n")
print(matches_m1)

cat("\nDiagnostics, M = 1:\n")
print(nn_match_diagnostics(matches_m1, n_total_treated = sum(dat$Z == 1)))

matches_m2 <- nn_match(
  data = dat,
  covariates = "X",
  M = 2,
  scaling = FALSE
)

cat("\nNearest-neighbor matches, M = 2:\n")
print(matches_m2)

cat("\nDiagnostics, M = 2:\n")
print(nn_match_diagnostics(matches_m2, n_total_treated = sum(dat$Z == 1)))
