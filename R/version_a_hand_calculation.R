# Fixed four-treated/four-control example for checking Version A by hand.
# The two CEM cells are X <= 0.5 and X > 0.5.

source(file.path("R", "few_treated_survival_methods.R"))

dat <- data.frame(
  id = c("t1", "t2", "t3", "t4", "c1", "c2", "c3", "c4"),
  Z = c(1, 1, 1, 1, 0, 0, 0, 0),
  X = c(0.10, 0.20, 1.10, 1.20, 0.15, 0.25, 1.05, 1.25),
  time = c(1, 2, 3, 4, 1, 2, 3, 5),
  status = c(1, 0, 1, 1, 1, 1, 0, 1)
)

cutpoints <- list(X = 0.5)
matched <- cem_keep_manual(
  dat, covariates = "X", cutpoints = cutpoints,
  treat_col = "Z", id_col = "id"
)

by <- baba_yoshida_cem_logrank(matched)
treated <- by$matched_data[by$matched_data$Z == 1, , drop = FALSE]
cell_groups <- stats::setNames(treated$cem_cell, as.character(treated$id))
test <- sign_change_test(
  S = by$contributions,
  groups = cell_groups[names(by$contributions)],
  alpha = 0.05,
  exact = TRUE,
  n_perm = 1000
)

cat("Version A hand-check example\n")
cat("CEM cells and matching status:\n")
print(matched[c("id", "Z", "X", "cem_cell", "cem_matched")], row.names = FALSE)
cat("\nEvent-time table (Y1, Y0, dN1, dN0, K, increment):\n")
print(by$event_table, row.names = FALSE)
cat("\nTreated-level contributions S_i:\n")
print(data.frame(id = names(by$contributions),
                 contribution = as.numeric(by$contributions)),
      row.names = FALSE)
cat("\nCell-level sign-change result:\n")
print(list(statistic = test$statistic, p_value = test$p_value,
           critical_value = test$critical_value, reject = test$reject,
           n_sign_units = test$n_sign_units,
           sign_grouping = "cem_cell"))

dir.create("results/00_current_results/05_version_a_hand_calculation",
           recursive = TRUE, showWarnings = FALSE)
write.csv(by$event_table,
          "results/00_current_results/05_version_a_hand_calculation/event_table.csv",
          row.names = FALSE)
write.csv(data.frame(
  id = names(by$contributions),
  contribution = as.numeric(by$contributions)
), "results/00_current_results/05_version_a_hand_calculation/treated_contributions.csv",
row.names = FALSE)
