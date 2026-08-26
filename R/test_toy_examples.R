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

expect_equal <- function(x, y, message, tolerance = 0) {
  if (tolerance == 0) {
    ok <- identical(x, y)
  } else {
    ok <- isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
  }
  if (!ok) {
    stop(message, call. = FALSE)
  }
}

quiet_value <- function(expr) {
  value <- NULL
  capture.output(value <- expr)
  value
}

run_toy_unit_tests <- function() {
  cat("Running toy/unit tests...\n")

  ## Basic helper functions
  ratio <- safe_ratio(c(1, 2, 3), c(1, 0, 2))
  expect_equal(ratio, c(1, 0, 1.5), "safe_ratio should return 0 when denominator is zero.")

  event_times <- event_times_until(
    time = c(3, 1, 2, 4, 2),
    status = c(1, 1, 0, 1, 1),
    tau = 3
  )
  expect_equal(event_times, c(1, 2, 3), "event_times_until returned the wrong sorted event times.")

  sc <- sign_change_test(c(1.5, -0.5, 0.8), alpha = 0.05)
  expect_true(is.finite(sc$statistic), "sign_change_test statistic is not finite.")
  expect_finite_p(sc$p_value, "sign_change_test")
  expect_true(sc$n_sign_units == 3, "sign_change_test reported the wrong number of sign units.")

  resolution_4 <- sign_change_resolution(4)
  resolution_5 <- sign_change_resolution(5)
  expect_true(resolution_4$critical_value_is_maximum, "Expected 4 sign units to have maximum critical value at alpha=0.05.")
  expect_true(resolution_5$smallest_positive_p_value == 1 / 32, "Unexpected p-value resolution for 5 sign units.")

  ## Matching helpers and diagnostics
  scalar_data <- toy_scalar_data()
  cov_matrix <- covariate_matrix(scalar_data, covariates = "X", scaling = FALSE)
  expect_true(nrow(cov_matrix) == nrow(scalar_data), "covariate_matrix row count changed.")
  expect_true(ncol(cov_matrix) == 1, "covariate_matrix should return one column for one covariate.")

  scalar_matches <- nn_match(scalar_data, covariates = "X", M = 1, scaling = FALSE)
  expect_true(all(c("treated_id", "control_id", "distance") %in% names(scalar_matches)), "nn_match output columns changed.")
  expect_true(nrow(scalar_matches) == 3, "NN toy example should return one control per treated.")

  scalar_diag <- nn_match_diagnostics(scalar_matches, n_total_treated = 3)
  expect_true(scalar_diag$n_matched_treated == 3, "NN toy example did not match all treated units.")
  expect_true(scalar_diag$n_discarded_treated == 0, "NN toy example should not discard treated units.")
  expect_true(scalar_diag$n_sign_units == 3, "NN toy example sign-unit diagnostic changed.")
  expect_true(scalar_diag$n_reused_controls >= 0, "NN reuse diagnostic should be nonnegative.")

  ferman_toy <- quiet_value(run_ferman_toy())
  expect_true(length(ferman_toy$contributions) == 3, "Ferman toy contribution count changed.")
  expect_finite_p(ferman_toy$p_value, "Ferman toy")

  ## CEM helpers and survival helpers
  survival_data <- toy_survival_data()
  cells <- cem_cells(survival_data, covariates = "cell_x", n_bins = 2)
  expect_true(length(cells) == nrow(survival_data), "cem_cells returned the wrong length.")
  expect_true(!any(is.na(cells)), "cem_cells should not return NA cells for the toy data.")

  cutpoints <- cem_internal_cutpoints(survival_data$X, n_bins = 2)
  expect_true(length(cutpoints) == 1, "cem_internal_cutpoints should return one cutpoint for two bins.")

  cem_data <- cem_keep(survival_data, covariates = "cell_x")
  cem_diag <- cem_match_diagnostics(cem_data)
  expect_true(cem_diag$n_total_treated == 3, "CEM toy total treated count changed.")
  expect_true(cem_diag$n_matched_treated == 3, "CEM toy example did not match all treated units.")
  expect_true(cem_diag$n_discarded_treated == 0, "CEM toy example should not discard treated units.")
  expect_true(cem_diag$n_sign_units == 2, "CEM toy sign-unit count changed.")

  weights <- cem_control_weights_at(cem_data, t = 3, treat_col = "Z", time_col = "time", cell_col = "cem_cell")
  expect_true(length(weights) == nrow(cem_data), "cem_control_weights_at should return one weight per subject.")
  expect_true(all(weights[cem_data$Z == 1] == 1), "cem_control_weights_at should give treated units weight one.")
  expect_true(all(is.finite(weights)), "cem_control_weights_at returned non-finite weights.")
  expect_true(all(weights >= 0), "cem_control_weights_at returned negative weights.")
  expect_true(weights[cem_data$id == "c6"] == 0, "cem_control_weights_at should give unmatched controls weight zero.")

  ## Method-level wrappers
  baba_toy <- quiet_value(run_cem_survival_toy())
  expect_true(is.finite(baba_toy$statistic), "Baba-Yoshida toy statistic is not finite.")
  expect_true(abs(sum(baba_toy$contributions) - baba_toy$statistic) < 1e-10, "Baba-Yoshida CEM contributions do not sum to statistic.")
  expect_finite_p(baba_toy$p_gaussian_two_sided, "Baba-Yoshida toy Gaussian")

  version_a <- quiet_value(run_version_a_toy())
  expect_true(version_a$version == "A_CEM_survival_sign_change", "Version A toy returned the wrong version label.")
  expect_finite_p(version_a$p_value, "Version A toy")
  expect_true(version_a$diagnostics$n_matched_treated == 3, "Version A toy matching diagnostic changed.")
  expect_true(version_a$n_contributions == version_a$diagnostics$n_matched_treated, "Version A contribution count should equal matched treated count.")
  expect_true(length(version_a$baba_yoshida$contributions) == version_a$n_contributions, "Version A nested contribution count changed.")
  expect_true(version_a$n_sign_units == version_a$diagnostics$n_sign_units, "Version A sign-unit diagnostic changed.")

  version_a_cell <- quiet_value(cem_survival_sign_test(
    data = toy_survival_data(),
    covariates = "cell_x",
    n_bins = 4,
    group_by_cell = TRUE,
    exact = TRUE,
    n_perm = 100
  ))
  expect_true(version_a_cell$version == "A_CEM_survival_sign_change", "Version A cell-level toy returned the wrong version label.")
  expect_finite_p(version_a_cell$p_value, "Version A cell-level toy")
  expect_true(version_a_cell$sign_grouping == "cem_cell", "Version A cell-level grouping was not used.")

  version_b <- quiet_value(run_version_b_toy())
  expect_true(version_b$version == "B_NN_survival_sign_change", "Version B toy returned the wrong version label.")
  expect_finite_p(version_b$p_value, "Version B toy")
  expect_true(version_b$diagnostics$n_matched_treated == 3, "Version B toy matching diagnostic changed.")
  expect_true(version_b$n_contributions == version_b$diagnostics$n_matched_treated, "Version B contribution count should equal matched treated count.")
  expect_true(length(version_b$nn_survival$contributions) == version_b$n_contributions, "Version B nested contribution count changed.")
  expect_true(version_b$n_sign_units == version_b$diagnostics$n_sign_units, "Version B sign-unit diagnostic changed.")

  ## Ferman alternative and tiny grid checks
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

  ferman_hd_dat <- generate_ferman_dgp(
    N1 = 3,
    N0 = 20,
    panel = "A",
    tau = 0,
    seed = 20260824,
    covariate_dim = 4
  )
  expect_true(
    all(paste0("X", 1:4) %in% names(ferman_hd_dat)),
    "Ferman high-dimensional DGP should create X1-X4."
  )
  ferman_hd <- run_one_ferman_replication(
    iter = 1,
    N1 = 3,
    N0 = 20,
    M = 1,
    panel = "A",
    tau = 0.5,
    scenario = "alternative",
    alternative_strength = "test",
    n_perm = 100,
    covariate_dim = 4
  )
  expect_true(ferman_hd$covariate_dim == 4, "Ferman high-dimensional replication did not record covariate_dim.")
  expect_finite_p(ferman_hd$p_value, "Ferman high-dimensional toy")

  ferman_grid <- run_ferman_replication_grid(
    n_iter = 1,
    panels = "A",
    N1_values = 5,
    N0 = 50,
    M_values = 1,
    covariate_dims = c(1L, 2L),
    tau_values = c(null = 0, alternative = 0.5),
    progress = FALSE,
    resume = FALSE
  )
  expect_true(nrow(ferman_grid$summary) == 4, "Ferman toy grid should have null/alternative rows for two covariate dimensions.")
  expect_true(all(c("null", "alternative") %in% ferman_grid$summary$scenario), "Ferman toy grid is missing null or alternative.")
  expect_true(all(c(1L, 2L) %in% ferman_grid$summary$covariate_dim), "Ferman toy grid is missing covariate dimensions.")

  cat("Toy/unit tests passed.\n")
  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_toy_unit_tests()
}
