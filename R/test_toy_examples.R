## Unit tests for toy examples and core method wrappers.
##
## Run from the project folder:
## Rscript R/test_toy_examples.R

find_project_root_for_tests <- function() {
  if (file.exists(file.path("R", "few_treated_survival_methods.R"))) {
    return(normalizePath("."))
  }
  if (file.exists("few_treated_survival_methods.R")) {
    return(normalizePath(".."))
  }
  stop("Run this script from the project folder or the R/ folder.")
}

project_root <- find_project_root_for_tests()
setwd(project_root)
source(file.path("R", "few_treated_survival_methods.R"))

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

expect_finite_p <- function(p, label) {
  expect_true(is.finite(p), paste(label, "p-value is not finite."))
  expect_true(p >= 0 && p <= 1, paste(label, "p-value is outside [0, 1]."))
}

quiet_value <- function(expr) {
  value <- NULL
  capture.output(value <- expr)
  value
}

run_toy_unit_tests <- function() {
  cat("Running toy/unit tests...\n")

  sc <- sign_change_test(c(1.5, -0.5, 0.8), alpha = 0.05)
  expect_true(is.finite(sc$statistic), "sign_change_test statistic is not finite.")
  expect_finite_p(sc$p_value, "sign_change_test")

  resolution_4 <- sign_change_resolution(4)
  resolution_5 <- sign_change_resolution(5)
  expect_true(resolution_4$critical_value_is_maximum, "Expected 4 sign units to have maximum critical value at alpha=0.05.")
  expect_true(resolution_5$smallest_positive_p_value == 1 / 32, "Unexpected p-value resolution for 5 sign units.")

  scalar_data <- toy_scalar_data()
  scalar_matches <- nn_match(scalar_data, covariates = "X", M = 1, scaling = FALSE)
  scalar_diag <- nn_match_diagnostics(scalar_matches, n_total_treated = 3)
  expect_true(scalar_diag$n_matched_treated == 3, "NN toy example did not match all treated units.")
  expect_true(scalar_diag$n_sign_units == 3, "NN toy example sign-unit diagnostic changed.")

  ferman_toy <- quiet_value(run_ferman_toy())
  expect_true(length(ferman_toy$contributions) == 3, "Ferman toy contribution count changed.")
  expect_finite_p(ferman_toy$p_value, "Ferman toy")

  survival_data <- toy_survival_data()
  cem_data <- cem_keep(survival_data, covariates = "cell_x")
  cem_diag <- cem_match_diagnostics(cem_data)
  expect_true(cem_diag$n_matched_treated == 3, "CEM toy example did not match all treated units.")

  baba_toy <- quiet_value(run_cem_survival_toy())
  expect_true(abs(sum(baba_toy$contributions) - baba_toy$statistic) < 1e-10, "Baba-Yoshida CEM contributions do not sum to statistic.")

  version_a <- quiet_value(run_version_a_toy())
  expect_true(version_a$version == "A_CEM_survival_sign_change", "Version A toy returned the wrong version label.")
  expect_finite_p(version_a$p_value, "Version A toy")

  version_b <- quiet_value(run_version_b_toy())
  expect_true(version_b$version == "B_NN_survival_sign_change", "Version B toy returned the wrong version label.")
  expect_finite_p(version_b$p_value, "Version B toy")

  ferman_alt <- run_one_ferman_replication(
    iter = 1,
    N1 = 5,
    N0 = 50,
    M = 1,
    panel = "A",
    tau = 0.5,
    scenario = "alternative",
    alternative_strength = "test",
    n_perm = 100
  )
  expect_true(ferman_alt$scenario == "alternative", "Ferman alternative scenario label was not recorded.")
  expect_true(abs(ferman_alt$tau - 0.5) < 1e-12, "Ferman alternative tau was not recorded.")
  expect_finite_p(ferman_alt$p_value, "Ferman alternative toy")

  ferman_grid <- run_ferman_replication_grid(
    n_iter = 1,
    panels = "A",
    N1_values = 5,
    N0 = 50,
    M_values = 1,
    tau_values = c(null = 0, alternative = 0.5),
    progress = FALSE,
    resume = FALSE
  )
  expect_true(nrow(ferman_grid$summary) == 2, "Ferman toy grid should have one null row and one alternative row.")
  expect_true(all(c("null", "alternative") %in% ferman_grid$summary$scenario), "Ferman toy grid is missing null or alternative.")

  cat("Toy/unit tests passed.\n")
  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_toy_unit_tests()
}
