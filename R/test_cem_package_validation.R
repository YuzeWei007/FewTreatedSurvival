## Compare manual CEM cells with package-based CEM from MatchIt.
##
## Run from either the project folder or the R/ folder:
## Rscript R/test_cem_package_validation.R

if (file.exists(file.path("R", "few_treated_survival_methods.R"))) {
  source(file.path("R", "few_treated_survival_methods.R"))
} else {
  source("few_treated_survival_methods.R")
}

if (!requireNamespace("MatchIt", quietly = TRUE)) {
  stop("Package 'MatchIt' is not installed.")
}

internal_cutpoints <- function(x, n_bins = 4) {
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  qs <- qs[is.finite(qs)]
  if (length(qs) <= 2) {
    return(numeric(0))
  }
  as.numeric(qs[-c(1, length(qs))])
}

compare_cem_engines <- function(data, covariates, n_bins = 4) {
  cutpoints <- stats::setNames(
    lapply(covariates, function(v) internal_cutpoints(data[[v]], n_bins = n_bins)),
    covariates
  )

  manual <- cem_keep(
    data = data,
    covariates = covariates,
    cutpoints = cutpoints,
    n_bins = n_bins,
    engine = "manual"
  )

  f <- stats::as.formula(paste("Z ~", paste(covariates, collapse = " + ")))
  pkg_fit <- MatchIt::matchit(
    formula = f,
    data = data,
    method = "cem",
    cutpoints = cutpoints
  )
  pkg <- MatchIt::match.data(pkg_fit, drop.unmatched = FALSE)
  pkg$matchit_matched <- pkg$weights > 0

  out <- data.frame(
    id = data$id,
    Z = data$Z,
    manual_cell = manual$cem_cell,
    manual_matched = manual$cem_matched,
    matchit_subclass = as.character(pkg$subclass),
    matchit_weight = pkg$weights,
    matchit_matched = pkg$matchit_matched,
    stringsAsFactors = FALSE
  )
  out$matched_agree <- out$manual_matched == out$matchit_matched

  list(
    cutpoints = cutpoints,
    comparison = out,
    all_matched_flags_agree = all(out$matched_agree)
  )
}

set.seed(20260811)
toy <- generate_survival_dgp(
  n1 = 10,
  n0 = 50,
  hr = 1,
  censor_rate = 0.08,
  x_shift_treated = 0
)

res <- compare_cem_engines(toy, covariates = "X", n_bins = 4)

cat("Cutpoints used by both engines:\n")
print(res$cutpoints)

cat("\nManual CEM vs MatchIt CEM matched flags:\n")
print(res$comparison)

cat("\nDo matched/unmatched flags agree? ", res$all_matched_flags_agree, "\n", sep = "")
