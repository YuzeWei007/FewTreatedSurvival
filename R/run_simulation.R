## Main command:
## Rscript R/run_simulation.R
##
## By default, this runs the complete 300-replication simulation:
## Baba-Yoshida survival N0=300 + Ferman N0=1000 + proposed Version A/B N0=1000.
##
## You can also override 300, e.g.:
## Rscript R/run_simulation.R complete_300 500

N_REPLICATIONS <- 300L
BABA_SURVIVAL_N0 <- 300L
BABA_SURVIVAL_N1_VALUES <- c(10L, 20L, 25L, 50L, 100L, 150L)
PROPOSED_SURVIVAL_N0 <- 1000L
PROPOSED_SURVIVAL_N1_VALUES <- c(5L, 10L, 25L, 50L)
SURVIVAL_M_VALUES <- c(1L, 2L)
SURVIVAL_OVERLAP_SHIFTS <- c(0, 0.75)
SURVIVAL_CENSOR_RATES <- c(0.04, 0.08)
SURVIVAL_ALTERNATIVE_HR_VALUES <- c(
  weak = 0.70,
  medium = 0.50,
  strong = 0.35
)
SURVIVAL_N_BINS <- function(n_total, n1, n0) {
  max(4L, as.integer(floor(n_total^(1 / 3))))
}
SURVIVAL_N_BINS_DESCRIPTION <- "max(4, floor((N0 + N1)^(1/3)))"
SURVIVAL_N_PERM <- 1000L
requested_cores <- suppressWarnings(as.integer(Sys.getenv("SIM_CORES", "2")))
if (is.na(requested_cores)) {
  requested_cores <- 2L
}
available_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(available_cores)) {
  available_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
}
if (is.na(available_cores)) {
  available_cores <- 2L
}
SURVIVAL_N_CORES <- min(
  max(1L, requested_cores),
  max(1L, available_cores)
)
FERMAN_N0 <- 1000L
FERMAN_N1_VALUES <- c(5L, 10L, 25L, 50L)
FERMAN_M_VALUES <- c(1L, 4L, 10L)
FERMAN_TAU_VALUES <- c(
  null = 0,
  alternative = 0.5
)

SURVIVAL_N0 <- BABA_SURVIVAL_N0
SURVIVAL_N1_VALUES <- BABA_SURVIVAL_N1_VALUES
SURVIVAL_N0_SENSITIVITY <- PROPOSED_SURVIVAL_N0
SURVIVAL_N1_SENSITIVITY_VALUES <- PROPOSED_SURVIVAL_N1_VALUES

find_project_root <- function() {
  if (file.exists(file.path("R", "few_treated_survival_methods.R"))) {
    return(normalizePath("."))
  }
  if (file.exists("few_treated_survival_methods.R")) {
    return(normalizePath(".."))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
    return(normalizePath(file.path(dirname(script_path), "..")))
  }

  stop("Cannot find project root. Run from the project folder or the R/ folder.")
}

project_root <- find_project_root()
setwd(project_root)
source(file.path("R", "few_treated_survival_methods.R"))

args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args) >= 1) args[[1]] else "complete_300"
n_iter <- if (length(args) >= 2) as.integer(args[[2]]) else N_REPLICATIONS

output_dir <- "results"
detail_dir <- file.path("results", "detailed")
grid_dir <- file.path("results", "grid_level")
run_timestamp <- Sys.getenv("SIM_RUN_ID", unset = "")
if (!nzchar(run_timestamp)) {
  run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
}
timestamped_run_dir <- file.path(output_dir, "runs", run_timestamp)

dir.create(output_dir, showWarnings = FALSE)
dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)

timestamped_result_path <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  if (startsWith(normalized, paste0(output_dir, "/"))) {
    relative <- sub(paste0("^", output_dir, "/"), "", normalized)
  } else {
    relative <- basename(normalized)
  }

  ext <- tools::file_ext(relative)
  stem <- if (nzchar(ext)) {
    sub(paste0("\\.", ext, "$"), "", relative)
  } else {
    relative
  }
  stamped_relative <- if (nzchar(ext)) {
    paste0(stem, "_", run_timestamp, ".", ext)
  } else {
    paste0(stem, "_", run_timestamp)
  }
  file.path(timestamped_run_dir, stamped_relative)
}

write_output_csv <- function(x, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, ...)

  stamped_file <- timestamped_result_path(file)
  dir.create(dirname(stamped_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, stamped_file, ...)
  invisible(stamped_file)
}

final_result_files <- c(
  "01_ferman_2021_nn_no_survival.csv",
  "02_baba_yoshida_2024_cem_survival.csv",
  "03_version_a_cem_survival_ferman_pvalue.csv",
  "04_version_b_nn_survival_ferman_pvalue.csv"
)

run_self_tests <- function() {
  quiet <- function(expr) {
    value <- NULL
    capture.output(value <- expr)
    value
  }

  sc <- sign_change_test(c(1.5, -0.5, 0.8), alpha = 0.05)
  stopifnot(is.finite(sc$statistic))
  stopifnot(sc$p_value >= 0 && sc$p_value <= 1)

  match_df <- data.frame(
    treated_id = c("t1", "t1", "t2", "t2", "t3"),
    control_id = c("c1", "c2", "c2", "c3", "c4")
  )
  cc <- connected_components_from_matches(match_df)
  stopifnot(cc[["t1"]] == cc[["t2"]])
  stopifnot(cc[["t3"]] != cc[["t1"]])

  dat_scalar <- toy_scalar_data()
  matches_scalar <- nn_match(dat_scalar, covariates = "X", M = 1, scaling = FALSE)
  diag_scalar <- nn_match_diagnostics(matches_scalar, n_total_treated = 3)
  stopifnot(diag_scalar$n_matched_treated == 3)
  stopifnot(diag_scalar$n_discarded_treated == 0)
  stopifnot(diag_scalar$n_sign_units == 3)

  dat_surv <- toy_survival_data()
  cem_data <- cem_keep(dat_surv, covariates = "cell_x")
  diag_cem <- cem_match_diagnostics(cem_data)
  stopifnot(diag_cem$n_total_treated == 3)
  stopifnot(diag_cem$n_matched_treated == 3)
  stopifnot(diag_cem$n_sign_units == 2)

  ferman <- quiet(run_ferman_toy())
  stopifnot(length(ferman$contributions) == 3)
  stopifnot(ferman$p_value >= 0 && ferman$p_value <= 1)

  by <- quiet(run_cem_survival_toy())
  stopifnot(abs(sum(by$contributions) - by$statistic) < 1e-10)

  version_a <- quiet(run_version_a_toy())
  stopifnot(version_a$version == "A_CEM_survival_sign_change")
  stopifnot(version_a$p_value >= 0 && version_a$p_value <= 1)

  version_b <- quiet(run_version_b_toy())
  stopifnot(version_b$version == "B_NN_survival_sign_change")
  stopifnot(version_b$p_value >= 0 && version_b$p_value <= 1)

  sim <- run_small_simulation(n_iter = 2, n1 = 5, n0 = 50, M = 1)
  required_cols <- c(
    "method", "rejection_rate", "mean_p_value", "p_value_q10",
    "p_value_q50", "p_value_q90", "mean_matched_treated",
    "mean_discarded_treated", "n_errors"
  )
  stopifnot(all(required_cols %in% names(sim$summary)))
  source(file.path("R", "test_toy_examples.R"))
  cat("All tests passed.\n")
}

write_small_simulation <- function() {
  dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)
  sim <- run_small_simulation(
    n_iter = ifelse(is.na(n_iter), 10L, n_iter),
    n1 = 5,
    n0 = 200,
    M = 2,
    alpha = 0.05,
    tau_scalar = 0,
    survival_hr = 1,
    n_bins = 4,
    seed_base = 20260726
  )

  write_output_csv(sim$results, file.path(detail_dir, "small_simulation_results.csv"), row.names = FALSE)
  write_output_csv(sim$summary, file.path(output_dir, "smoke_simulation_summary.csv"), row.names = FALSE)
  print(sim$summary)
}

write_method_grid_files <- function(summary, prefix) {
  write_output_csv(
    summary,
    file.path(grid_dir, paste0(prefix, "_grid_summary_all_methods.csv")),
    row.names = FALSE
  )

  for (method_name in sort(unique(summary$method))) {
    method_rows <- summary[summary$method == method_name, ]
    safe_method_name <- gsub("[^A-Za-z0-9]+", "_", method_name)
    write_output_csv(
      method_rows,
      file.path(grid_dir, paste0(prefix, "_", safe_method_name, "_grid_summary.csv")),
      row.names = FALSE
    )
  }

  method_summary <- aggregate(
    cbind(rejection_rate, mean_p_value, mean_matched_treated, n_errors) ~
      scenario + alternative_strength + survival_hr + method,
    data = summary,
    FUN = mean,
    na.rm = TRUE
  )
  method_summary$settings_count <- aggregate(
    n_rows ~ scenario + alternative_strength + survival_hr + method,
    data = summary,
    FUN = length
  )$n_rows
  method_summary <- method_summary[
    order(method_summary$scenario, method_summary$alternative_strength, method_summary$method),
  ]
  write_output_csv(
    method_summary,
    file.path(grid_dir, paste0(prefix, "_method_comparison_summary.csv")),
    row.names = FALSE
  )

  n1_summary <- aggregate(
    cbind(rejection_rate, mean_p_value, mean_matched_treated, n_errors) ~
      scenario + alternative_strength + survival_hr + n1 + n_bins + method,
    data = summary,
    FUN = mean,
    na.rm = TRUE
  )
  n1_summary <- n1_summary[
    order(n1_summary$scenario, n1_summary$alternative_strength, n1_summary$method, n1_summary$n1),
  ]
  write_output_csv(
    n1_summary,
    file.path(grid_dir, paste0(prefix, "_by_n1_summary.csv")),
    row.names = FALSE
  )

  strength_n1_summary <- aggregate(
    cbind(rejection_rate, mean_p_value, mean_matched_treated, n_errors) ~
      scenario + alternative_strength + survival_hr + n1 + n_bins + method,
    data = summary,
    FUN = mean,
    na.rm = TRUE
  )
  strength_n1_summary <- strength_n1_summary[
    order(
      strength_n1_summary$scenario,
      strength_n1_summary$alternative_strength,
      strength_n1_summary$method,
      strength_n1_summary$n1
    ),
  ]
  write_output_csv(
    strength_n1_summary,
    file.path(grid_dir, paste0(prefix, "_by_strength_n1_summary.csv")),
    row.names = FALSE
  )
}

dedupe_ferman_summary <- function(summary) {
  key_cols <- c("scenario", "alternative_strength", "tau", "panel", "N1", "M")
  if (!all(key_cols %in% names(summary))) {
    return(summary)
  }
  summary <- summary[order(
    summary$scenario,
    summary$alternative_strength,
    summary$tau,
    summary$panel,
    summary$N1,
    summary$M
  ), ]
  summary[!duplicated(summary[key_cols]), , drop = FALSE]
}

write_ferman_table2_null_rows <- function(summary, prefix, alpha_value) {
  null_rows <- summary[summary$scenario == "null", , drop = FALSE]
  null_rows <- null_rows[order(null_rows$panel, null_rows$N1, null_rows$M), ]
  null_rows$paper_table_reference <- "Ferman 2021 Table 2"
  null_rows$comparison_note <- paste0(
    "Table-2-aligned run: alpha=", alpha_value,
    "; compare null rejection_rate by panel, N1, and M."
  )
  write_output_csv(
    null_rows,
    file.path(grid_dir, paste0(prefix, "_table2_null_rows.csv")),
    row.names = FALSE
  )
  invisible(null_rows)
}

write_formal_grid <- function(default_iter = 20L,
                              n0_value = SURVIVAL_N0,
                              n1_values = SURVIVAL_N1_VALUES,
                              methods = c(
                                "baba_yoshida_cem_gaussian",
                                "version_a_cem_sign",
                                "version_b_nn_sign"
                              ),
                              prefix = "survival_n0_300",
                              include_scalar_ferman = FALSE,
                              update_main_outputs = TRUE,
                              n_cores = 1L) {
  dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)
  methods <- match.arg(
    methods,
    choices = c(
      "baba_yoshida_cem_gaussian",
      "version_a_cem_sign",
      "version_b_nn_sign"
    ),
    several.ok = TRUE
  )
  actual_iter <- ifelse(is.na(n_iter), default_iter, n_iter)
  cat("Running survival formal grid with", actual_iter, "replications per grid cell.\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n")
  cat("Replication cores:", n_cores, "\n")
  cat("Survival grid: N0=", n0_value,
      ", N1={", paste(n1_values, collapse = ", "), "}",
      ", M={", paste(SURVIVAL_M_VALUES, collapse = ", "), "}",
      ", overlap_shift={", paste(SURVIVAL_OVERLAP_SHIFTS, collapse = ", "), "}",
      ", censor_rate={", paste(SURVIVAL_CENSOR_RATES, collapse = ", "), "}",
      ", CEM n_bins=", SURVIVAL_N_BINS_DESCRIPTION,
      ", alternative HR={", paste(names(SURVIVAL_ALTERNATIVE_HR_VALUES), SURVIVAL_ALTERNATIVE_HR_VALUES, sep = "=", collapse = ", "), "}",
      ", sign-change n_perm=", SURVIVAL_N_PERM, ".\n",
      sep = "")
  grid <- run_formal_simulation_grid(
    n_iter = actual_iter,
    n1_values = n1_values,
    n0 = n0_value,
    M_values = SURVIVAL_M_VALUES,
    overlap_shifts = SURVIVAL_OVERLAP_SHIFTS,
    censor_rates = SURVIVAL_CENSOR_RATES,
    scenarios = c("null", "alternative"),
    alternative_hr_values = SURVIVAL_ALTERNATIVE_HR_VALUES,
    alpha = 0.05,
    n_bins = SURVIVAL_N_BINS,
    seed_base = 20260728,
    progress = TRUE,
    include_scalar_ferman = include_scalar_ferman,
    n_perm = SURVIVAL_N_PERM,
    methods = methods,
    checkpoint_dir = grid_dir,
    checkpoint_prefix = prefix,
    resume = TRUE,
    n_cores = n_cores
  )

  if (update_main_outputs) {
    write_output_csv(grid$results, file.path(detail_dir, "formal_simulation_grid_results.csv"), row.names = FALSE)
    write_output_csv(grid$summary, file.path(output_dir, "formal_simulation_grid_summary.csv"), row.names = FALSE)
  }
  write_output_csv(grid$results, file.path(detail_dir, paste0(prefix, "_raw_results.csv")), row.names = FALSE)
  write_output_csv(grid$summary, file.path(grid_dir, paste0(prefix, "_grid_summary_all_methods.csv")), row.names = FALSE)
  write_method_grid_files(grid$summary, prefix)
  print(grid$summary)
  if (update_main_outputs) {
    cat("Main formal grid summary written to results/formal_simulation_grid_summary.csv\n")
  }
  cat("Full method-level survival grid files written under results/grid_level/ with prefix ", prefix, ".\n", sep = "")
}

write_complete_simulation <- function(default_iter = 300L) {
  actual_iter <- ifelse(is.na(n_iter), default_iter, n_iter)
  cat("Running COMPLETE simulation with", actual_iter, "replications per grid cell.\n")
  cat("Step 1/3: Baba-Yoshida CEM survival baseline, N0=300.\n")
  write_formal_grid(
    default_iter = actual_iter,
    n0_value = BABA_SURVIVAL_N0,
    n1_values = BABA_SURVIVAL_N1_VALUES,
    methods = "baba_yoshida_cem_gaussian",
    prefix = "baba_yoshida_strength_n0_300",
    update_main_outputs = TRUE
  )

  cat("Step 2/3: Ferman non-survival replication, N0=1000, null plus alternative.\n")
  write_ferman_replication(default_iter = actual_iter, include_alternative = TRUE)

  cat("Step 3/3: proposed Version A/B survival methods, N0=1000.\n")
  write_formal_grid(
    default_iter = actual_iter,
    n0_value = PROPOSED_SURVIVAL_N0,
    n1_values = PROPOSED_SURVIVAL_N1_VALUES,
    methods = c("version_a_cem_sign", "version_b_nn_sign"),
    prefix = "version_ab_strength_n0_1000",
    update_main_outputs = FALSE
  )

  write_report_tables()
}

write_version_ab_30min_pilot <- function(default_iter = 20L) {
  actual_iter <- if (length(args) >= 2) as.integer(args[[2]]) else default_iter
  pilot_n_perm <- 300L
  prefix <- "version_ab_30min_pilot_n0_1000"
  dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Running 30-minute Version A/B pilot.\n")
  cat("Purpose: quick trend check, not the full 300-replication grid.\n")
  cat("Replications per grid cell:", actual_iter, "\n")
  cat("Methods: version_a_cem_sign, version_b_nn_sign\n")
  cat("N0: 1000\n")
  cat("N1: 5, 10, 25, 50\n")
  cat("Scenarios: null and strong alternative only\n")
  cat("M: 1\n")
  cat("overlap_shift: 0\n")
  cat("censor_rate: 0.08\n")
  cat("Sign-change n_perm:", pilot_n_perm, "\n")
  cat("Replication cores:", SURVIVAL_N_CORES, "\n")

  grid <- run_formal_simulation_grid(
    n_iter = actual_iter,
    n1_values = c(5L, 10L, 25L, 50L),
    n0 = PROPOSED_SURVIVAL_N0,
    M_values = c(1L),
    overlap_shifts = c(0),
    censor_rates = c(0.08),
    scenarios = c("null", "alternative"),
    alternative_hr_values = c(strong = SURVIVAL_ALTERNATIVE_HR_VALUES[["strong"]]),
    alpha = 0.05,
    n_bins = SURVIVAL_N_BINS,
    seed_base = 20260813,
    progress = TRUE,
    include_scalar_ferman = FALSE,
    n_perm = pilot_n_perm,
    methods = c("version_a_cem_sign", "version_b_nn_sign"),
    checkpoint_dir = grid_dir,
    checkpoint_prefix = prefix,
    resume = TRUE,
    n_cores = SURVIVAL_N_CORES
  )

  write_output_csv(grid$results, file.path(detail_dir, paste0(prefix, "_raw_results.csv")), row.names = FALSE)
  write_output_csv(grid$summary, file.path(grid_dir, paste0(prefix, "_grid_summary_all_methods.csv")), row.names = FALSE)
  write_method_grid_files(grid$summary, prefix)
  print(grid$summary)
  cat("30-minute pilot files written under results/grid_level/ with prefix ", prefix, ".\n", sep = "")
}

write_ferman_replication <- function(default_iter = 5L,
                                     include_alternative = TRUE,
                                     alpha_value = 0.05,
                                     prefix_suffix = NULL,
                                     write_table2_null = FALSE) {
  dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)
  actual_iter <- ifelse(is.na(n_iter), default_iter, n_iter)
  tau_values <- if (isTRUE(include_alternative)) FERMAN_TAU_VALUES else c(null = 0)
  base_prefix <- if (isTRUE(include_alternative)) {
    "ferman_2021_nn_no_survival_with_alternative_120_settings"
  } else {
    "ferman_2021_nn_no_survival_60_settings"
  }
  prefix <- if (is.null(prefix_suffix)) base_prefix else paste0(base_prefix, "_", prefix_suffix)
  cat("Running Ferman replication grid with", actual_iter, "replications per grid cell.\n")
  cat("Ferman scenarios:", paste(names(tau_values), collapse = ", "), "\n")
  cat("Ferman alpha:", alpha_value, "\n")
  pilot <- run_ferman_replication_grid(
    n_iter = actual_iter,
    panels = c("A", "B", "C", "D", "E"),
    N1_values = FERMAN_N1_VALUES,
    M_values = FERMAN_M_VALUES,
    N0 = FERMAN_N0,
    tau_values = tau_values,
    alpha = alpha_value,
    n_perm = 1000,
    seed_base = 20260728,
    progress = TRUE,
    checkpoint_dir = grid_dir,
    checkpoint_prefix = prefix,
    resume = TRUE
  )
  pilot$summary <- dedupe_ferman_summary(pilot$summary)

  write_output_csv(pilot$results, file.path(detail_dir, "ferman_replication_results.csv"), row.names = FALSE)
  write_output_csv(pilot$summary, file.path(output_dir, "ferman_replication_summary.csv"), row.names = FALSE)
  write_output_csv(pilot$summary, file.path(grid_dir, paste0(prefix, ".csv")), row.names = FALSE)
  write_output_csv(pilot$summary, file.path(grid_dir, paste0(prefix, "_grid_summary.csv")), row.names = FALSE)
  if (isTRUE(write_table2_null)) {
    write_ferman_table2_null_rows(pilot$summary, prefix, alpha_value)
  }

  if (file.exists(file.path(output_dir, "old_attempts_ferman_saved_results_summary.csv"))) {
    null_summary <- pilot$summary[pilot$summary$scenario == "null", , drop = FALSE]
    comparison <- compare_to_old_ferman_outputs(null_summary)
    write_output_csv(comparison, file.path(output_dir, "ferman_replication_vs_old_attempts.csv"), row.names = FALSE)
    write_output_csv(comparison, file.path(grid_dir, paste0(prefix, "_null_vs_old_attempts.csv")), row.names = FALSE)
    print(pilot$summary)
    cat("\nNull-only comparison to old saved outputs:\n")
    print(comparison)
  } else {
    print(pilot$summary)
  }
}

print_key_results <- function() {
  method_path <- file.path(output_dir, "method_comparison_summary.csv")
  censor_path <- file.path(output_dir, "censoring_sensitivity_summary.csv")
  ferman_path <- file.path(output_dir, "ferman_panel_comparison_summary.csv")

  if (file.exists(method_path)) {
    cat("\n=== Survival method comparison ===\n")
    print(read.csv(method_path, stringsAsFactors = FALSE))
  }

  if (file.exists(censor_path)) {
    cat("\n=== Survival censoring sensitivity ===\n")
    print(read.csv(censor_path, stringsAsFactors = FALSE))
  }

  if (file.exists(ferman_path)) {
    cat("\n=== Ferman panel comparison ===\n")
    print(read.csv(ferman_path, stringsAsFactors = FALSE))
  }
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

markdown_table <- function(df) {
  if (nrow(df) == 0) {
    return("_No rows available._")
  }

  df[] <- lapply(df, as.character)
  header <- paste(names(df), collapse = " | ")
  divider <- paste(rep("---", ncol(df)), collapse = " | ")
  rows <- apply(df, 1, paste, collapse = " | ")
  paste(c(header, divider, rows), collapse = "\n")
}

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

print_run_plan <- function() {
  n_survival_scenarios <- 1 + length(SURVIVAL_ALTERNATIVE_HR_VALUES)
  baba_settings <- n_survival_scenarios *
    length(BABA_SURVIVAL_N1_VALUES) *
    length(SURVIVAL_M_VALUES) *
    length(SURVIVAL_OVERLAP_SHIFTS) *
    length(SURVIVAL_CENSOR_RATES)
  proposed_settings <- n_survival_scenarios *
    length(PROPOSED_SURVIVAL_N1_VALUES) *
    length(SURVIVAL_M_VALUES) *
    length(SURVIVAL_OVERLAP_SHIFTS) *
    length(SURVIVAL_CENSOR_RATES)
  proposed_methods <- c("version_a_cem_sign", "version_b_nn_sign")
  ferman_settings <- length(FERMAN_TAU_VALUES) *
    length(c("A", "B", "C", "D", "E")) *
    length(FERMAN_N1_VALUES) *
    length(FERMAN_M_VALUES)

  cat("Complete run plan\n")
  cat("=================\n")
  cat("Default task: complete_300\n")
  cat("Default replications per grid cell:", N_REPLICATIONS, "\n\n")

  cat("1. Baba-Yoshida CEM survival baseline\n")
  cat("   Methods: baba_yoshida_cem_gaussian\n")
  cat("   N0:", BABA_SURVIVAL_N0, "\n")
  cat("   N1:", paste(BABA_SURVIVAL_N1_VALUES, collapse = ", "), "\n")
  cat("   M:", paste(SURVIVAL_M_VALUES, collapse = ", "), "\n")
  cat("   overlap_shift:", paste(SURVIVAL_OVERLAP_SHIFTS, collapse = ", "), "\n")
  cat("   censor_rate:", paste(SURVIVAL_CENSOR_RATES, collapse = ", "), "\n")
  cat("   CEM n_bins:", SURVIVAL_N_BINS_DESCRIPTION, "\n")
  cat("   alternative HR:", paste(names(SURVIVAL_ALTERNATIVE_HR_VALUES), SURVIVAL_ALTERNATIVE_HR_VALUES, sep = "=", collapse = ", "), "\n")
  cat("   grid rows:", baba_settings, "\n")
  cat("   output prefix: results/grid_level/baba_yoshida_strength_n0_300_*\n\n")

  cat("2. Ferman non-survival replication\n")
  cat("   N0:", FERMAN_N0, "\n")
  cat("   N1:", paste(FERMAN_N1_VALUES, collapse = ", "), "\n")
  cat("   M:", paste(FERMAN_M_VALUES, collapse = ", "), "\n")
  cat("   panels: A, B, C, D, E\n")
  cat("   scenarios:", paste(names(FERMAN_TAU_VALUES), collapse = ", "), "\n")
  cat("   alternative tau:", FERMAN_TAU_VALUES[["alternative"]], "\n")
  cat("   grid rows:", ferman_settings, "\n")
  cat("   output: results/grid_level/ferman_2021_nn_no_survival_with_alternative_120_settings.csv\n\n")

  cat("3. Proposed survival methods\n")
  cat("   Methods:", paste(proposed_methods, collapse = ", "), "\n")
  cat("   N0:", PROPOSED_SURVIVAL_N0, "\n")
  cat("   N1:", paste(PROPOSED_SURVIVAL_N1_VALUES, collapse = ", "), "\n")
  cat("   CEM n_bins:", SURVIVAL_N_BINS_DESCRIPTION, "\n")
  cat("   alternative HR:", paste(names(SURVIVAL_ALTERNATIVE_HR_VALUES), SURVIVAL_ALTERNATIVE_HR_VALUES, sep = "=", collapse = ", "), "\n")
  cat("   grid settings per method:", proposed_settings, "\n")
  cat("   all-method grid rows:", proposed_settings * length(proposed_methods), "\n")
  cat("   output prefix: results/grid_level/version_ab_strength_n0_1000_*\n\n")

  cat("Run command:\n")
  cat("   Rscript R/run_simulation.R complete_300\n")
}

write_final_result_files <- function() {
  method_summary <- read_if_exists(file.path(output_dir, "method_comparison_summary.csv"))
  ferman_detail <- read_if_exists(file.path(output_dir, "ferman_replication_vs_old_attempts.csv"))
  ferman_summary <- read_if_exists(file.path(output_dir, "ferman_replication_summary.csv"))
  ferman_panel <- read_if_exists(file.path(output_dir, "ferman_panel_comparison_summary.csv"))
  baba_grid <- read_if_exists(file.path(
    grid_dir,
    "baba_yoshida_strength_n0_300_baba_yoshida_cem_gaussian_grid_summary.csv"
  ))
  version_a_grid <- read_if_exists(file.path(
    grid_dir,
    "version_ab_strength_n0_1000_version_a_cem_sign_grid_summary.csv"
  ))
  version_b_grid <- read_if_exists(file.path(
    grid_dir,
    "version_ab_strength_n0_1000_version_b_nn_sign_grid_summary.csv"
  ))

  if (!is.null(ferman_detail)) {
    write_output_csv(
      ferman_detail,
      file.path(output_dir, "01_ferman_2021_nn_no_survival.csv"),
      row.names = FALSE
    )
  } else if (!is.null(ferman_summary)) {
    write_output_csv(
      ferman_summary,
      file.path(output_dir, "01_ferman_2021_nn_no_survival.csv"),
      row.names = FALSE
    )
  } else if (!is.null(ferman_panel)) {
    write_output_csv(
      ferman_panel,
      file.path(output_dir, "01_ferman_2021_nn_no_survival.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(baba_grid)) {
    write_output_csv(
      baba_grid,
      file.path(output_dir, "02_baba_yoshida_2024_cem_survival.csv"),
      row.names = FALSE
    )
  } else if (!is.null(method_summary)) {
    rows <- method_summary[method_summary$method == "baba_yoshida_cem_gaussian", ]
    if (nrow(rows) > 0) {
      write_output_csv(
        rows,
        file.path(output_dir, "02_baba_yoshida_2024_cem_survival.csv"),
        row.names = FALSE
      )
    }
  }

  if (!is.null(version_a_grid)) {
    write_output_csv(
      version_a_grid,
      file.path(output_dir, "03_version_a_cem_survival_ferman_pvalue.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(version_b_grid)) {
    write_output_csv(
      version_b_grid,
      file.path(output_dir, "04_version_b_nn_survival_ferman_pvalue.csv"),
      row.names = FALSE
    )
  }
}

clean_result_folder <- function() {
  keep <- file.path(output_dir, final_result_files)
  files <- list.files(output_dir, full.names = TRUE, recursive = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir == FALSE]
  extra_files <- setdiff(normalizePath(files, mustWork = FALSE), normalizePath(keep, mustWork = FALSE))

  if (length(extra_files) > 0) {
    unlink(extra_files)
  }
  if (dir.exists(detail_dir)) {
    unlink(detail_dir, recursive = TRUE)
  }
}

write_report_tables <- function() {
  small_summary <- read_if_exists(file.path(output_dir, "smoke_simulation_summary.csv"))
  formal_summary <- read_if_exists(file.path(output_dir, "formal_simulation_grid_summary.csv"))
  ferman_summary <- read_if_exists(file.path(output_dir, "ferman_replication_summary.csv"))
  ferman_compare <- read_if_exists(file.path(output_dir, "ferman_replication_vs_old_attempts.csv"))
  if (is.null(ferman_compare)) {
    ferman_compare <- read_if_exists(file.path(output_dir, "ferman_replication_pilot_vs_old_attempts.csv"))
  }

  formal_method_summary <- NULL
  formal_censor_summary <- NULL
  if (!is.null(formal_summary)) {
    formal_method_summary <- aggregate(
      cbind(rejection_rate, mean_p_value, mean_matched_treated, n_errors) ~
        scenario + alternative_strength + survival_hr + method,
      data = formal_summary,
      FUN = mean,
      na.rm = TRUE
    )
    formal_method_summary$settings_count <- aggregate(
      n_rows ~ scenario + alternative_strength + survival_hr + method,
      data = formal_summary,
      FUN = length
    )$n_rows
    formal_method_summary <- formal_method_summary[
      order(formal_method_summary$scenario, formal_method_summary$alternative_strength, formal_method_summary$method),
    ]
    write_output_csv(
      formal_method_summary,
      file.path(output_dir, "method_comparison_summary.csv"),
      row.names = FALSE
    )

    formal_censor_summary <- aggregate(
      cbind(rejection_rate, mean_p_value, mean_matched_treated, n_errors) ~
        scenario + alternative_strength + survival_hr + censor_rate + method,
      data = formal_summary,
      FUN = mean,
      na.rm = TRUE
    )
    formal_censor_summary <- formal_censor_summary[
      order(
        formal_censor_summary$scenario,
        formal_censor_summary$alternative_strength,
        formal_censor_summary$censor_rate,
        formal_censor_summary$method
      ),
    ]
    write_output_csv(
      formal_censor_summary,
      file.path(output_dir, "censoring_sensitivity_summary.csv"),
      row.names = FALSE
    )
  }

  ferman_panel_summary <- NULL
  if (!is.null(ferman_compare)) {
    ferman_compare$abs_rejection_rate_diff <- abs(ferman_compare$rejection_rate_diff)
    ferman_panel_summary <- aggregate(
      cbind(rejection_rate, old_attempt_rejection_rate, abs_rejection_rate_diff,
            mean_sign_units, mean_reused_controls) ~ panel,
      data = ferman_compare,
      FUN = mean,
      na.rm = TRUE
    )
    names(ferman_panel_summary)[names(ferman_panel_summary) == "abs_rejection_rate_diff"] <- "mean_abs_diff_vs_old_attempt"
    ferman_panel_summary <- ferman_panel_summary[order(ferman_panel_summary$panel), ]
    write_output_csv(
      ferman_panel_summary,
      file.path(output_dir, "ferman_panel_comparison_summary.csv"),
      row.names = FALSE
    )
  } else if (!is.null(ferman_summary)) {
    ferman_panel_summary <- aggregate(
      cbind(rejection_rate, mean_sign_units, mean_reused_controls) ~ panel,
      data = ferman_summary,
      FUN = mean,
      na.rm = TRUE
    )
    ferman_panel_summary <- ferman_panel_summary[order(ferman_panel_summary$panel), ]
    write_output_csv(
      ferman_panel_summary,
      file.path(output_dir, "ferman_panel_comparison_summary.csv"),
      row.names = FALSE
    )
  }

  lines <- c(
    "# Simulation Results Summary",
    "",
    "This report is generated from the current simulation outputs. The `n_iter` columns in the CSV files record the number of replications used.",
    "",
    "## Generated CSV Files",
    "",
    "- `results/method_comparison_summary.csv`: formal survival grid collapsed by scenario, alternative strength, hazard ratio, and method.",
    "- `results/censoring_sensitivity_summary.csv`: formal survival grid collapsed by scenario, alternative strength, hazard ratio, censoring rate, and method.",
    "- `results/ferman_panel_comparison_summary.csv`: Ferman replication compared with old saved GitHub outputs, collapsed by panel.",
    "- `results/grid_level/`: full grid-level summaries. These files are kept for checking every method and every simulation setting.",
    "",
    "## Small Simulation",
    "",
    "The smoke simulation checks that the four method wrappers all return p-values, rejection indicators, and matching diagnostics."
  )

  if (!is.null(small_summary)) {
    small_display <- small_summary[, intersect(
      c("method", "rejection_rate", "mean_p_value", "mean_sign_units",
        "mean_matched_treated", "mean_discarded_treated", "n_errors"),
      names(small_summary)
    )]
    num_cols <- vapply(small_display, is.numeric, logical(1))
    small_display[num_cols] <- lapply(small_display[num_cols], fmt_num)
    lines <- c(lines, "", markdown_table(small_display))
  }

  lines <- c(
    lines,
    "",
    "## Formal Survival Grid",
    "",
    "The formal survival grid reports null plus weak, medium, and strong alternatives separately. The alternative is controlled by the survival hazard ratio: weak HR=0.70, medium HR=0.50, and strong HR=0.35. CEM uses the dynamic bin rule `max(4, floor((N0 + N1)^(1/3)))`. Baba-Yoshida uses `N1 = 10, 20, 25, 50, 100, 150`; Version A/B use `N1 = 5, 10, 25, 50`. The non-survival Ferman replication is reported separately."
  )

  if (!is.null(formal_method_summary)) {
    formal_display <- formal_method_summary
    names(formal_display)[names(formal_display) == "rejection_rate"] <- "avg_rejection_rate"
    names(formal_display)[names(formal_display) == "mean_p_value"] <- "avg_mean_p_value"
    names(formal_display)[names(formal_display) == "mean_matched_treated"] <- "avg_matched_treated"
    names(formal_display)[names(formal_display) == "n_errors"] <- "avg_errors"
    num_cols <- vapply(formal_display, is.numeric, logical(1))
    formal_display[num_cols] <- lapply(formal_display[num_cols], fmt_num)
    lines <- c(lines, "", markdown_table(formal_display))
  }

  lines <- c(
    lines,
    "",
    "## Ferman Replication",
    "",
    "The Ferman runner is a low-dependency reimplementation of the relevant nearest-neighbor/sign-change logic. It compares rejection rates against old saved GitHub outputs when those comparison targets are present."
  )

  if (!is.null(ferman_panel_summary)) {
    ferman_display <- ferman_panel_summary
    if ("old_attempt_rejection_rate" %in% names(ferman_display)) {
      names(ferman_display) <- c(
        "panel", "pilot_rejection_rate", "old_attempt_rejection_rate",
        "mean_abs_diff_vs_old_attempt", "mean_sign_units", "mean_reused_controls"
      )
    }
    num_cols <- vapply(ferman_display, is.numeric, logical(1))
    ferman_display[num_cols] <- lapply(ferman_display[num_cols], fmt_num)
    lines <- c(lines, "", markdown_table(ferman_display))
  }

  lines <- c(
    lines,
    "",
    "## Interpretation",
    "",
    "- Completed: the two matching methods and the two p-value construction paths are implemented as runnable modules.",
    "- Completed: Version A and Version B are both included in the pilot simulation grid.",
    "- Completed: Baba-Yoshida Gaussian CEM inference uses the paper-style Proposition 3.2 Week 1 variance estimate in the reported Gaussian p-value.",
    "- Completed: Ferman replication has a standalone pilot runner and comparison target extracted from old saved GitHub outputs.",
    "- Remaining caveat: exact paper-level reproduction still depends on final DGP/grid confirmation and whether exact CSMatch dependency installation is required."
  )

  writeLines(lines, "results/simulation_report.md")
  write_final_result_files()
  cat("Wrote final result files and kept full grid-level outputs in results/grid_level/.\n")
  print_key_results()
}

if (task == "help") {
  cat("Usage: Rscript R/run_simulation.R <task> [n_iter]\n")
  cat("Tasks: plan, tests, unit_tests, small, formal, ferman, ferman_null, survival_300, survival_n0_1000_300, survival_n0_1000_300_parallel, version_ab_30min_pilot, ferman_300, ferman_table2_300, ferman_null_300, complete_300, paper_300, survival_final, survival_n0_1000_final, ferman_final, final, report, all\n")
} else if (task == "plan") {
  print_run_plan()
} else if (task == "tests") {
  run_self_tests()
} else if (task == "unit_tests") {
  source(file.path("R", "test_toy_examples.R"))
} else if (task == "small") {
  write_small_simulation()
} else if (task == "formal") {
  write_formal_grid(default_iter = 20L)
} else if (task == "ferman") {
  write_ferman_replication(default_iter = 5L, include_alternative = TRUE)
} else if (task == "ferman_null") {
  write_ferman_replication(default_iter = 5L, include_alternative = FALSE)
} else if (task == "survival_300") {
  write_formal_grid(
    default_iter = 300L,
    n0_value = BABA_SURVIVAL_N0,
    n1_values = BABA_SURVIVAL_N1_VALUES,
    methods = "baba_yoshida_cem_gaussian",
    prefix = "baba_yoshida_strength_n0_300"
  )
  write_report_tables()
} else if (task == "survival_n0_1000_300") {
  write_formal_grid(
    default_iter = 300L,
    n0_value = PROPOSED_SURVIVAL_N0,
    n1_values = PROPOSED_SURVIVAL_N1_VALUES,
    methods = c("version_a_cem_sign", "version_b_nn_sign"),
    prefix = "version_ab_strength_n0_1000",
    update_main_outputs = FALSE
  )
} else if (task == "survival_n0_1000_300_parallel") {
  write_formal_grid(
    default_iter = 300L,
    n0_value = PROPOSED_SURVIVAL_N0,
    n1_values = PROPOSED_SURVIVAL_N1_VALUES,
    methods = c("version_a_cem_sign", "version_b_nn_sign"),
    prefix = "version_ab_strength_n0_1000",
    update_main_outputs = FALSE,
    n_cores = SURVIVAL_N_CORES
  )
} else if (task == "version_ab_30min_pilot") {
  write_version_ab_30min_pilot(default_iter = 20L)
} else if (task == "ferman_300") {
  write_ferman_replication(default_iter = 300L, include_alternative = TRUE)
  write_report_tables()
} else if (task == "ferman_table2_300") {
  write_ferman_replication(
    default_iter = 300L,
    include_alternative = TRUE,
    alpha_value = 0.10,
    prefix_suffix = "table2_alpha10",
    write_table2_null = TRUE
  )
} else if (task == "ferman_null_300") {
  write_ferman_replication(default_iter = 300L, include_alternative = FALSE)
  write_report_tables()
} else if (task == "complete_300") {
  write_complete_simulation(default_iter = 300L)
} else if (task == "paper_300") {
  write_complete_simulation(default_iter = 300L)
} else if (task == "survival_final") {
  write_formal_grid(
    default_iter = 300L,
    n0_value = BABA_SURVIVAL_N0,
    n1_values = BABA_SURVIVAL_N1_VALUES,
    methods = "baba_yoshida_cem_gaussian",
    prefix = "baba_yoshida_strength_n0_300"
  )
} else if (task == "survival_n0_1000_final") {
  write_formal_grid(
    default_iter = 300L,
    n0_value = PROPOSED_SURVIVAL_N0,
    n1_values = PROPOSED_SURVIVAL_N1_VALUES,
    methods = c("version_a_cem_sign", "version_b_nn_sign"),
    prefix = "version_ab_strength_n0_1000",
    update_main_outputs = FALSE
  )
} else if (task == "ferman_final") {
  write_ferman_replication(default_iter = 300L, include_alternative = TRUE)
} else if (task == "report") {
  write_report_tables()
} else if (task == "final") {
  write_complete_simulation(default_iter = 300L)
} else if (task == "all") {
  run_self_tests()
  write_small_simulation()
  write_formal_grid(default_iter = 20L)
  write_ferman_replication(default_iter = 5L, include_alternative = TRUE)
  write_report_tables()
} else {
  stop("Unknown task: ", task)
}
