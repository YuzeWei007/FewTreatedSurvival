# Few Treated Survival Project - Main Code
#
# This single source file is organized into clear modules:
# 1. p-value construction
# 2. matching and diagnostics
# 3. survival contributions and tests
# 4. toy examples
# 5. simulation grids
# 6. Ferman replication

# ============================================================
# MODULE 1: P-value construction master
#
# This block contains the Ferman-style sign-change/randomization
# p-value machinery used by the scalar Ferman replication and by
# proposed survival Version A/B.
# ============================================================

safe_ratio <- function(num, den) {
  ifelse(den == 0, 0, num / den)
}

treated_t_stat <- function(S) {
  S <- as.numeric(S)
  n <- length(S)
  if (n < 2) {
    stop("Need at least two treated-level contributions.")
  }
  centered_ss <- sum((S - mean(S))^2)
  if (centered_ss == 0) {
    return(ifelse(mean(S) == 0, 0, Inf))
  }
  abs(mean(S)) / sqrt(centered_ss / (n - 1))
}

.sign_group_cache <- new.env(parent = emptyenv())

all_group_signs <- function(n_groups) {
  key <- as.character(n_groups)
  if (exists(key, envir = .sign_group_cache, inherits = FALSE)) {
    return(get(key, envir = .sign_group_cache, inherits = FALSE))
  }
  grid <- expand.grid(rep(list(c(-1, 1)), n_groups))
  signs <- as.matrix(grid)
  assign(key, signs, envir = .sign_group_cache)
  signs
}

treated_t_stat_many <- function(signed_values) {
  n <- ncol(signed_values)
  means <- rowMeans(signed_values)
  centered <- signed_values - means
  centered_ss <- rowSums(centered * centered)
  out <- abs(means) / sqrt(centered_ss / (n - 1))
  zero_var <- centered_ss == 0
  if (any(zero_var)) {
    out[zero_var] <- ifelse(means[zero_var] == 0, 0, Inf)
  }
  out
}

sign_matrix_from_groups <- function(groups, exact = TRUE, n_perm = 10000,
                                    exact_max_groups = 20, seed = NULL) {
  groups <- as.factor(groups)
  group_levels <- levels(groups)
  n_groups <- length(group_levels)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (exact && n_groups <= exact_max_groups) {
    group_signs <- all_group_signs(n_groups)
  } else {
    group_signs <- matrix(
      sample(c(-1, 1), n_perm * n_groups, replace = TRUE),
      nrow = n_perm,
      ncol = n_groups
    )
  }

  signs <- group_signs[, as.integer(groups), drop = FALSE]
  colnames(signs) <- names(groups)
  signs
}

sign_change_test <- function(S, groups = NULL, alpha = 0.05, exact = TRUE,
                             n_perm = 10000, exact_max_groups = 20,
                             seed = NULL) {
  S <- as.numeric(S)
  if (is.null(groups)) {
    groups <- seq_along(S)
  }
  if (length(groups) != length(S)) {
    stop("groups must be NULL or have the same length as S.")
  }

  T_obs <- treated_t_stat(S)
  signs <- sign_matrix_from_groups(
    groups = groups,
    exact = exact,
    n_perm = n_perm,
    exact_max_groups = exact_max_groups,
    seed = seed
  )
  T_null <- treated_t_stat_many(sweep(signs, 2, S, "*"))
  T_sorted <- sort(T_null)
  k <- ceiling((1 - alpha) * length(T_sorted))
  critical_value <- T_sorted[k]
  p_value <- mean(T_null >= T_obs)

  list(
    statistic = T_obs,
    p_value = p_value,
    critical_value = critical_value,
    reject = T_obs > critical_value,
    null_statistics = T_null,
    n_sign_units = length(unique(groups)),
    n_contributions = length(S),
    alpha = alpha,
    exact = exact && length(unique(groups)) <= exact_max_groups
  )
}

sign_change_resolution <- function(n_sign_units, alpha = 0.05) {
  if (length(n_sign_units) != 1 || is.na(n_sign_units) || n_sign_units < 1) {
    stop("n_sign_units must be one positive integer.")
  }
  n_sign_units <- as.integer(n_sign_units)
  n_assignments <- 2^n_sign_units
  k <- ceiling((1 - alpha) * n_assignments)
  data.frame(
    n_sign_units = n_sign_units,
    n_assignments = n_assignments,
    alpha = alpha,
    critical_index = k,
    smallest_positive_p_value = 1 / n_assignments,
    critical_value_is_maximum = k >= n_assignments,
    stringsAsFactors = FALSE
  )
}

connected_components_from_matches <- function(matched_controls) {
  if (is.data.frame(matched_controls)) {
    split_controls <- split(matched_controls$control_id, matched_controls$treated_id)
  } else if (is.matrix(matched_controls)) {
    split_controls <- lapply(seq_len(nrow(matched_controls)), function(i) matched_controls[i, ])
    names(split_controls) <- rownames(matched_controls)
  } else if (is.list(matched_controls)) {
    split_controls <- matched_controls
  } else {
    stop("matched_controls must be a data.frame, matrix, or list.")
  }

  ids <- names(split_controls)
  n <- length(split_controls)
  parent <- seq_len(n)

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union_roots <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) {
      parent[rb] <<- ra
    }
  }

  clean_controls <- lapply(split_controls, function(x) unique(stats::na.omit(as.character(x))))
  if (n >= 2) {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        if (length(intersect(clean_controls[[i]], clean_controls[[j]])) > 0) {
          union_roots(i, j)
        }
      }
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  groups <- match(roots, unique(roots))
  names(groups) <- ids
  groups
}


# ============================================================
# MODULE 2: Matching and matching-diagnostics master
#
# This block contains nearest-neighbor matching, CEM matching,
# package/manual CEM wrappers, and reuse/sign-unit diagnostics.
# ============================================================

covariate_matrix <- function(data, covariates, scaling = TRUE) {
  X <- as.matrix(data[, covariates, drop = FALSE])
  storage.mode(X) <- "double"

  if (isTRUE(scaling)) {
    sds <- apply(X, 2, stats::sd)
    sds[sds == 0 | is.na(sds)] <- 1
    X <- sweep(X, 2, sds, "/")
  } else if (is.numeric(scaling) && length(scaling) > 1) {
    X <- sweep(X, 2, scaling, "*")
  } else if (is.numeric(scaling) && length(scaling) == 1 && scaling != 1) {
    X <- X * scaling
  }

  X
}

nn_match <- function(data, covariates, treat_col = "Z", id_col = "id",
                     M = 1, scaling = TRUE) {
  if (!(id_col %in% names(data))) {
    data[[id_col]] <- as.character(seq_len(nrow(data)))
  }

  treated_idx <- which(data[[treat_col]] == 1 | data[[treat_col]] == TRUE)
  control_idx <- which(!(data[[treat_col]] == 1 | data[[treat_col]] == TRUE))
  if (length(treated_idx) == 0 || length(control_idx) == 0) {
    stop("Need at least one treated and one control unit.")
  }
  if (M > length(control_idx)) {
    stop("M cannot exceed the number of controls.")
  }

  X <- covariate_matrix(data, covariates, scaling = scaling)
  Xt <- X[treated_idx, , drop = FALSE]
  Xc <- X[control_idx, , drop = FALSE]

  rows <- lapply(seq_along(treated_idx), function(i) {
    diffs <- sweep(Xc, 2, Xt[i, ], "-")
    distances <- sqrt(rowSums(diffs^2))
    ord <- order(distances)[seq_len(M)]
    data.frame(
      treated_id = as.character(data[[id_col]][treated_idx[i]]),
      control_id = as.character(data[[id_col]][control_idx[ord]]),
      rank = seq_len(M),
      distance = distances[ord],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

coarsen_one <- function(x, cutpoints = NULL, n_bins = NULL) {
  if (!is.numeric(x)) {
    return(as.character(x))
  }

  if (is.null(cutpoints)) {
    if (is.null(n_bins)) {
      n_bins <- 4
    }
    probs <- seq(0, 1, length.out = n_bins + 1)
    cutpoints <- unique(stats::quantile(x, probs = probs, na.rm = TRUE))
  } else {
    cutpoints <- sort(unique(c(-Inf, cutpoints, Inf)))
  }

  as.character(cut(x, breaks = cutpoints, include.lowest = TRUE, right = TRUE))
}

cem_cells <- function(data, covariates, cutpoints = NULL, n_bins = 4) {
  pieces <- lapply(covariates, function(v) {
    cp <- if (!is.null(cutpoints) && v %in% names(cutpoints)) cutpoints[[v]] else NULL
    coarsen_one(data[[v]], cutpoints = cp, n_bins = n_bins)
  })
  do.call(paste, c(pieces, sep = " | "))
}

cem_internal_cutpoints <- function(x, n_bins = 4) {
  if (!is.numeric(x)) {
    return(NULL)
  }
  qs <- unique(stats::quantile(
    x,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  ))
  qs <- qs[is.finite(qs)]
  if (length(qs) <= 2) {
    return(numeric(0))
  }
  as.numeric(qs[-c(1, length(qs))])
}

cem_matchit_cutpoints <- function(data, covariates, cutpoints = NULL, n_bins = 4) {
  stats::setNames(lapply(covariates, function(v) {
    if (!is.null(cutpoints) && v %in% names(cutpoints)) {
      cp <- as.numeric(cutpoints[[v]])
      cp <- cp[is.finite(cp)]
      return(sort(unique(cp)))
    }
    cem_internal_cutpoints(data[[v]], n_bins = n_bins)
  }), covariates)
}

cem_keep_manual <- function(data, covariates, treat_col = "Z", id_col = "id",
                            cutpoints = NULL, n_bins = 4) {
  if (!(id_col %in% names(data))) {
    data[[id_col]] <- as.character(seq_len(nrow(data)))
  }

  data$cem_cell <- cem_cells(data, covariates, cutpoints = cutpoints, n_bins = n_bins)
  is_treated <- data[[treat_col]] == 1 | data[[treat_col]] == TRUE
  tab <- table(data$cem_cell, is_treated)
  keep_cells <- rownames(tab)[tab[, "TRUE"] > 0 & tab[, "FALSE"] > 0]
  data$cem_matched <- data$cem_cell %in% keep_cells
  data$cem_engine <- "manual"
  data$cem_weight <- ifelse(data$cem_matched, 1, 0)
  data
}

cem_keep_matchit <- function(data, covariates, treat_col = "Z", id_col = "id",
                             cutpoints = NULL, n_bins = 4) {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is not installed.")
  }
  if (!all(vapply(data[, covariates, drop = FALSE], is.numeric, logical(1)))) {
    stop("MatchIt CEM path currently expects numeric covariates.")
  }
  if (!(id_col %in% names(data))) {
    data[[id_col]] <- as.character(seq_len(nrow(data)))
  }

  row_col <- ".cem_matchit_row_id"
  while (row_col %in% names(data)) {
    row_col <- paste0(row_col, "_")
  }
  match_data <- data
  match_data[[row_col]] <- seq_len(nrow(match_data))

  cp <- cem_matchit_cutpoints(
    data = match_data,
    covariates = covariates,
    cutpoints = cutpoints,
    n_bins = n_bins
  )
  f <- stats::as.formula(paste(treat_col, "~", paste(covariates, collapse = " + ")))
  fit <- MatchIt::matchit(
    formula = f,
    data = match_data,
    method = "cem",
    cutpoints = cp
  )
  matched <- MatchIt::match.data(fit, drop.unmatched = FALSE)
  matched <- matched[order(matched[[row_col]]), , drop = FALSE]

  out <- data
  out$cem_cell <- cem_cells(out, covariates, cutpoints = cp, n_bins = n_bins)
  if ("subclass" %in% names(matched)) {
    subclass <- as.character(matched$subclass)
    out$cem_cell <- ifelse(is.na(subclass), out$cem_cell, paste0("matchit_", subclass))
  }
  out$cem_matched <- matched$weights > 0
  out$cem_engine <- "MatchIt"
  out$cem_weight <- matched$weights
  out
}

cem_keep <- function(data, covariates, treat_col = "Z", id_col = "id",
                     cutpoints = NULL, n_bins = 4,
                     engine = c("auto", "matchit", "manual")) {
  engine <- match.arg(engine)
  if (engine == "manual") {
    return(cem_keep_manual(
      data = data,
      covariates = covariates,
      treat_col = treat_col,
      id_col = id_col,
      cutpoints = cutpoints,
      n_bins = n_bins
    ))
  }

  if (engine == "auto") {
    can_use_matchit <- requireNamespace("MatchIt", quietly = TRUE) &&
      all(vapply(data[, covariates, drop = FALSE], is.numeric, logical(1)))
    if (!can_use_matchit) {
      return(cem_keep_manual(
        data = data,
        covariates = covariates,
        treat_col = treat_col,
        id_col = id_col,
        cutpoints = cutpoints,
        n_bins = n_bins
      ))
    }
  }

  cem_keep_matchit(
    data = data,
    covariates = covariates,
    treat_col = treat_col,
    id_col = id_col,
    cutpoints = cutpoints,
    n_bins = n_bins
  )
}

nn_match_diagnostics <- function(matches, n_total_treated = NULL) {
  split_matches <- split(matches$control_id, matches$treated_id)
  groups <- connected_components_from_matches(matches)
  control_counts <- table(matches$control_id)
  control_reuse <- as.integer(control_counts)

  n_matched_treated <- length(split_matches)
  if (is.null(n_total_treated)) {
    n_total_treated <- n_matched_treated
  }

  data.frame(
    n_total_treated = n_total_treated,
    n_matched_treated = n_matched_treated,
    n_discarded_treated = n_total_treated - n_matched_treated,
    n_sign_units = length(unique(groups)),
    n_controls_used = length(control_counts),
    number_unique_controls = length(control_counts),
    total_control_uses = nrow(matches),
    mean_controls_per_treated = mean(lengths(split_matches)),
    min_control_reuse = min(control_reuse),
    mean_control_reuse = mean(control_reuse),
    max_control_reuse = max(control_reuse),
    n_reused_controls = sum(control_reuse > 1),
    unmatched_treated = n_total_treated - n_matched_treated,
    unmatched_controls = NA_integer_,
    reuse_definition = "nn_matched_index_count_per_control",
    stringsAsFactors = FALSE
  )
}

cem_match_diagnostics <- function(cem_data, treat_col = "Z", cell_col = "cem_cell") {
  is_treated <- cem_data[[treat_col]] == 1 | cem_data[[treat_col]] == TRUE
  matched <- cem_data$cem_matched %in% TRUE
  treated_cells <- unique(cem_data[[cell_col]][matched & is_treated])
  matched_controls <- matched & !is_treated
  control_cells <- cem_data[[cell_col]][matched_controls]
  treated_per_cell <- tapply(as.numeric(is_treated & matched), cem_data[[cell_col]], sum)
  baseline_reuse <- as.numeric(treated_per_cell[control_cells])
  baseline_reuse <- baseline_reuse[!is.na(baseline_reuse) & baseline_reuse > 0]
  n_matched_treated <- sum(is_treated & matched)
  n_matched_controls <- sum(matched_controls)

  data.frame(
    n_total_treated = sum(is_treated),
    n_matched_treated = n_matched_treated,
    n_discarded_treated = sum(is_treated & !matched),
    unmatched_treated = sum(is_treated & !matched),
    n_sign_units = length(treated_cells),
    n_matched_controls = n_matched_controls,
    n_controls_used = n_matched_controls,
    number_unique_controls = n_matched_controls,
    unmatched_controls = sum(!is_treated & !matched),
    n_matched_cells = length(unique(cem_data[[cell_col]][matched])),
    min_control_reuse = ifelse(length(baseline_reuse) == 0, NA_real_, min(baseline_reuse)),
    mean_control_reuse = ifelse(length(baseline_reuse) == 0, NA_real_, mean(baseline_reuse)),
    max_control_reuse = ifelse(length(baseline_reuse) == 0, NA_real_, max(baseline_reuse)),
    n_reused_controls = ifelse(length(baseline_reuse) == 0, NA_integer_, sum(baseline_reuse > 1)),
    reuse_definition = "cem_baseline_cell_treated_count_per_matched_control",
    stringsAsFactors = FALSE
  )
}

ferman_scalar_contributions <- function(data, matches, outcome_col = "Y",
                                        id_col = "id", tau0 = 0) {
  ids <- as.character(data[[id_col]])
  names(ids) <- seq_along(ids)
  y_by_id <- stats::setNames(data[[outcome_col]], as.character(data[[id_col]]))
  split_matches <- split(matches$control_id, matches$treated_id)

  S <- vapply(names(split_matches), function(tid) {
    y_by_id[[tid]] - mean(y_by_id[as.character(split_matches[[tid]])]) - tau0
  }, numeric(1))

  S
}

ferman_scalar_nn_test <- function(data, covariates, outcome_col = "Y",
                                  treat_col = "Z", id_col = "id", M = 1,
                                  tau0 = 0, alpha = 0.05, scaling = TRUE,
                                  exact = TRUE, n_perm = 10000, seed = NULL) {
  matches <- nn_match(
    data = data,
    covariates = covariates,
    treat_col = treat_col,
    id_col = id_col,
    M = M,
    scaling = scaling
  )
  S <- ferman_scalar_contributions(
    data = data,
    matches = matches,
    outcome_col = outcome_col,
    id_col = id_col,
    tau0 = tau0
  )
  groups <- connected_components_from_matches(matches)
  test <- sign_change_test(
    S = S,
    groups = groups[names(S)],
    alpha = alpha,
    exact = exact,
    n_perm = n_perm,
    seed = seed
  )
  test$matches <- matches
  test$contributions <- S
  n_total_treated <- sum(data[[treat_col]] == 1 | data[[treat_col]] == TRUE)
  test$diagnostics <- nn_match_diagnostics(matches, n_total_treated = n_total_treated)
  test
}


# ============================================================
# MODULE 3: Survival contribution master
#
# This block turns matched survival data into treated-level
# contributions for Baba-Yoshida CEM, Version A, and Version B.
# ============================================================

event_times_until <- function(time, status, tau = Inf) {
  sort(unique(time[status == 1 & time <= tau]))
}

cem_control_weights_at <- function(data, t, treat_col, time_col, cell_col) {
  is_treated <- data[[treat_col]] == 1 | data[[treat_col]] == TRUE
  at_risk <- data[[time_col]] >= t
  cells <- data[[cell_col]]

  treated_risk <- tapply(as.numeric(is_treated & at_risk), cells, sum)
  control_risk <- tapply(as.numeric(!is_treated & at_risk), cells, sum)
  cell_ratio <- safe_ratio(treated_risk[names(control_risk)], control_risk)
  cell_ratio[is.na(cell_ratio)] <- 0

  weights <- numeric(nrow(data))
  control_rows <- which(!is_treated)
  weights[is_treated] <- 1
  weights[control_rows] <- cell_ratio[cells[control_rows]]
  weights[is.na(weights)] <- 0
  weights
}

baba_yoshida_cem_logrank <- function(data, time_col = "time",
                                     status_col = "status",
                                     treat_col = "Z",
                                     cell_col = "cem_cell",
                                     id_col = "id",
                                     tau = Inf) {
  data <- data[data$cem_matched %in% TRUE, , drop = FALSE]
  if (nrow(data) == 0) {
    stop("No CEM-matched observations.")
  }
  if (!(id_col %in% names(data))) {
    data[[id_col]] <- as.character(seq_len(nrow(data)))
  }

  is_treated <- data[[treat_col]] == 1 | data[[treat_col]] == TRUE
  treated_ids <- as.character(data[[id_col]][is_treated])
  event_times <- event_times_until(data[[time_col]], data[[status_col]], tau)

  w0 <- cem_control_weights_at(data, 0, treat_col, time_col, cell_col)
  Y10 <- sum(w0 * is_treated * (data[[time_col]] >= 0))
  Y00 <- sum(w0 * (!is_treated) * (data[[time_col]] >= 0))

  contributions <- stats::setNames(rep(0, length(treated_ids)), treated_ids)
  W <- 0
  event_table <- data.frame()

  for (tk in event_times) {
    weights <- cem_control_weights_at(data, tk, treat_col, time_col, cell_col)
    at_risk <- data[[time_col]] >= tk
    events <- data[[time_col]] == tk & data[[status_col]] == 1

    Y1 <- sum(weights * is_treated * at_risk)
    Y0 <- sum(weights * (!is_treated) * at_risk)
    dN1 <- sum(weights * is_treated * events)
    dN0 <- sum(weights * (!is_treated) * events)

    K <- sqrt(safe_ratio(Y10 + Y00, Y10 * Y00)) * safe_ratio(Y1 * Y0, Y1 + Y0)
    increment <- K * (safe_ratio(dN1, Y1) - safe_ratio(dN0, Y0))
    W <- W + increment

    treated_rows <- which(is_treated)
    for (r in treated_rows) {
      tid <- as.character(data[[id_col]][r])
      dNi <- as.numeric(events[r])
      Yi <- as.numeric(at_risk[r])
      contributions[tid] <- contributions[tid] +
        K * (safe_ratio(dNi, Y1) - safe_ratio(Yi, Y1) * safe_ratio(dN0, Y0))
    }

    event_table <- rbind(
      event_table,
      data.frame(time = tk, Y1 = Y1, Y0 = Y0, dN1 = dN1, dN0 = dN0,
                 K = K, increment = increment)
    )
  }

  treated_event_count <- sum(is_treated & data[[status_col]] == 1 & data[[time_col]] <= tau)
  variance_week1_simple <- treated_event_count / length(treated_ids)^2
  variance_paper_prop32_w1 <- 0.5 * treated_event_count / length(treated_ids)
  Z <- safe_ratio(W, sqrt(variance_paper_prop32_w1))

  list(
    statistic = W,
    contributions = contributions,
    event_table = event_table,
    variance_simple = variance_week1_simple,
    variance_week1_simple = variance_week1_simple,
    variance_paper_prop32_w1 = variance_paper_prop32_w1,
    z_gaussian = Z,
    p_gaussian_two_sided = 2 * (1 - stats::pnorm(abs(Z))),
    matched_data = data
  )
}

cem_survival_sign_test <- function(data, covariates, time_col = "time",
                                   status_col = "status", treat_col = "Z",
                                   id_col = "id", cutpoints = NULL,
                                   n_bins = 4, tau = Inf, alpha = 0.05,
                                   group_by_cell = FALSE, exact = TRUE,
                                   n_perm = 10000, seed = NULL) {
  cem_data <- cem_keep(
    data = data,
    covariates = covariates,
    treat_col = treat_col,
    id_col = id_col,
    cutpoints = cutpoints,
    n_bins = n_bins
  )
  diagnostics <- cem_match_diagnostics(cem_data, treat_col = treat_col, cell_col = "cem_cell")
  by <- baba_yoshida_cem_logrank(
    data = cem_data,
    time_col = time_col,
    status_col = status_col,
    treat_col = treat_col,
    cell_col = "cem_cell",
    id_col = id_col,
    tau = tau
  )

  treated_data <- by$matched_data[by$matched_data[[treat_col]] == 1 |
                                    by$matched_data[[treat_col]] == TRUE, ]
  groups <- if (group_by_cell) treated_data$cem_cell else treated_data[[id_col]]
  names(groups) <- as.character(treated_data[[id_col]])

  test <- sign_change_test(
    S = by$contributions,
    groups = groups[names(by$contributions)],
    alpha = alpha,
    exact = exact,
    n_perm = n_perm,
    seed = seed
  )

  test$version <- "A_CEM_survival_sign_change"
  test$sign_grouping <- if (group_by_cell) "cem_cell" else "treated_patient"
  test$baba_yoshida <- by
  test$cem_data <- cem_data
  test$diagnostics <- diagnostics
  test
}

nn_survival_contributions <- function(data, matches, time_col = "time",
                                      status_col = "status", treat_col = "Z",
                                      id_col = "id", tau = Inf,
                                      q_fun = function(t) 1) {
  if (!(id_col %in% names(data))) {
    data[[id_col]] <- as.character(seq_len(nrow(data)))
  }
  row_by_id <- stats::setNames(seq_len(nrow(data)), as.character(data[[id_col]]))
  split_matches <- split(matches$control_id, matches$treated_id)
  event_times <- event_times_until(data[[time_col]], data[[status_col]], tau)

  S <- stats::setNames(rep(0, length(split_matches)), names(split_matches))
  detail <- data.frame()

  for (tid in names(split_matches)) {
    member_ids <- c(tid, as.character(split_matches[[tid]]))
    member_rows <- row_by_id[member_ids]
    treated_row <- row_by_id[[tid]]

    for (tk in event_times) {
      at_risk <- data[[time_col]][member_rows] >= tk
      events <- data[[time_col]][member_rows] == tk & data[[status_col]][member_rows] == 1
      Ri <- sum(at_risk)
      Di <- sum(events)
      Yi <- as.numeric(data[[time_col]][treated_row] >= tk)
      dNi <- as.numeric(data[[time_col]][treated_row] == tk &&
                          data[[status_col]][treated_row] == 1)
      increment <- q_fun(tk) * (dNi - safe_ratio(Yi * Di, Ri))
      S[tid] <- S[tid] + increment

      detail <- rbind(
        detail,
        data.frame(treated_id = tid, time = tk, Ri = Ri, Di = Di,
                   Yi = Yi, dNi = dNi, increment = increment)
      )
    }
  }

  list(contributions = S, detail = detail)
}

nn_survival_sign_test <- function(data, covariates, time_col = "time",
                                  status_col = "status", treat_col = "Z",
                                  id_col = "id", M = 1, scaling = TRUE,
                                  tau = Inf, alpha = 0.05, exact = TRUE,
                                  n_perm = 10000, seed = NULL) {
  matches <- nn_match(
    data = data,
    covariates = covariates,
    treat_col = treat_col,
    id_col = id_col,
    M = M,
    scaling = scaling
  )
  nn <- nn_survival_contributions(
    data = data,
    matches = matches,
    time_col = time_col,
    status_col = status_col,
    treat_col = treat_col,
    id_col = id_col,
    tau = tau
  )
  groups <- connected_components_from_matches(matches)
  test <- sign_change_test(
    S = nn$contributions,
    groups = groups[names(nn$contributions)],
    alpha = alpha,
    exact = exact,
    n_perm = n_perm,
    seed = seed
  )
  test$version <- "B_NN_survival_sign_change"
  test$sign_grouping <- "shared_control_connected_component"
  test$matches <- matches
  test$nn_survival <- nn
  n_total_treated <- sum(data[[treat_col]] == 1 | data[[treat_col]] == TRUE)
  test$diagnostics <- nn_match_diagnostics(matches, n_total_treated = n_total_treated)
  test
}


# ============================================================
# MODULE 4: Toy examples master
#
# This block gives tiny scalar and survival datasets that are used
# for quick checks and unit tests.
# ============================================================

toy_scalar_data <- function() {
  data.frame(
    id = c("t1", "t2", "t3", "c1", "c2", "c3", "c4", "c5"),
    Z = c(1, 1, 1, 0, 0, 0, 0, 0),
    X = c(0.0, 1.0, 2.0, 0.1, 0.9, 1.8, 2.2, 3.0),
    Y = c(4.0, 5.0, 8.0, 2.0, 4.0, 5.0, 7.0, 9.0)
  )
}

toy_survival_data <- function() {
  data.frame(
    id = c("t1", "t2", "t3", "c1", "c2", "c3", "c4", "c5", "c6"),
    Z = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
    X = c(0.2, 0.4, 1.2, 0.1, 0.3, 0.5, 1.0, 1.3, 2.0),
    cell_x = c("low", "low", "high", "low", "low", "low", "high", "high", "far"),
    time = c(5, 7, 4, 6, 3, 8, 5, 9, 2),
    status = c(1, 0, 1, 1, 1, 0, 0, 1, 1)
  )
}

run_ferman_toy <- function() {
  dat <- toy_scalar_data()
  res <- ferman_scalar_nn_test(
    data = dat,
    covariates = "X",
    outcome_col = "Y",
    M = 1,
    alpha = 0.05,
    scaling = FALSE
  )
  print(res[c("statistic", "p_value", "critical_value", "reject", "n_sign_units")])
  print(res$contributions)
  invisible(res)
}

run_cem_survival_toy <- function() {
  dat <- toy_survival_data()
  dat$cem_cell <- dat$cell_x
  dat$cem_matched <- dat$cell_x %in% c("low", "high")
  by <- baba_yoshida_cem_logrank(dat)
  print(by$event_table)
  print(by$contributions)
  print(by[c("statistic", "variance_week1_simple", "variance_paper_prop32_w1",
             "z_gaussian", "p_gaussian_two_sided")])
  invisible(by)
}

run_version_a_toy <- function() {
  dat <- toy_survival_data()
  res <- cem_survival_sign_test(
    data = dat,
    covariates = "cell_x",
    n_bins = 2,
    group_by_cell = TRUE,
    alpha = 0.05
  )
  print(res[c("version", "statistic", "p_value", "critical_value",
              "reject", "n_sign_units")])
  print(res$baba_yoshida$contributions)
  invisible(res)
}

run_version_b_toy <- function() {
  dat <- toy_survival_data()
  res <- nn_survival_sign_test(
    data = dat,
    covariates = "X",
    M = 2,
    scaling = FALSE,
    alpha = 0.05
  )
  print(res[c("version", "statistic", "p_value", "critical_value",
              "reject", "n_sign_units")])
  print(res$nn_survival$contributions)
  invisible(res)
}

run_all_toy_examples <- function() {
  list(
    ferman_scalar = run_ferman_toy(),
    baba_yoshida_cem = run_cem_survival_toy(),
    version_a = run_version_a_toy(),
    version_b = run_version_b_toy()
  )
}


# ============================================================
# MODULE 5: Simulation grid master
#
# This block generates scalar/survival DGPs, runs one replication,
# summarizes results, and manages checkpointed survival grids.
# ============================================================

generate_scalar_dgp <- function(n1 = 5, n0 = 200, tau = 0, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- n1 + n0
  Z <- c(rep(1, n1), rep(0, n0))
  X <- c(stats::rnorm(n1, 0, 1), stats::rnorm(n0, 0, 1))
  mu0 <- X
  Y0 <- mu0 + stats::rnorm(n)
  Y1 <- Y0 + tau
  Y <- ifelse(Z == 1, Y1, Y0)

  data.frame(
    id = paste0(ifelse(Z == 1, "t", "c"), seq_len(n)),
    Z = Z,
    X = X,
    Y = Y,
    Y0 = Y0,
    Y1 = Y1
  )
}

generate_survival_dgp <- function(n1 = 5, n0 = 200, hr = 1,
                                  censor_rate = 0.08, x_shift_treated = 0,
                                  seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- n1 + n0
  Z <- c(rep(1, n1), rep(0, n0))
  X <- c(stats::rnorm(n1, x_shift_treated, 1), stats::rnorm(n0, 0, 1))
  base_rate <- 0.08 * exp(0.35 * X)
  event_rate <- base_rate * ifelse(Z == 1, hr, 1)
  event_time <- stats::rexp(n, rate = event_rate)
  censor_time <- stats::rexp(n, rate = censor_rate)
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)

  data.frame(
    id = paste0(ifelse(Z == 1, "t", "c"), seq_len(n)),
    Z = Z,
    X = X,
    time = time,
    status = status
  )
}

safe_method_result <- function(expr, method) {
  tryCatch(
    expr,
    error = function(e) {
      list(method = method, error = e$message)
    }
  )
}

run_one_simulation <- function(iter = 1, n1 = 5, n0 = 200, M = 2,
                               alpha = 0.05, tau_scalar = 0,
                               survival_hr = 1, n_bins = 4,
                               censor_rate = 0.08, x_shift_treated = 0,
                               seed_base = 1000,
                               include_scalar_ferman = TRUE,
                               n_perm = 1000,
                               methods = c(
                                 "baba_yoshida_cem_gaussian",
                                 "version_a_cem_sign",
                                 "version_b_nn_sign"
                               )) {
  methods <- match.arg(
    methods,
    choices = c(
      "baba_yoshida_cem_gaussian",
      "version_a_cem_sign",
      "version_b_nn_sign"
    ),
    several.ok = TRUE
  )
  scalar <- NULL
  if (include_scalar_ferman) {
    scalar <- generate_scalar_dgp(
      n1 = n1,
      n0 = n0,
      tau = tau_scalar,
      seed = seed_base + iter
    )
  }
  survival <- generate_survival_dgp(
    n1 = n1,
    n0 = n0,
    hr = survival_hr,
    censor_rate = censor_rate,
    x_shift_treated = x_shift_treated,
    seed = seed_base + 100000 + iter
  )

  ferman <- NULL
  if (include_scalar_ferman) {
    ferman <- safe_method_result(
      ferman_scalar_nn_test(
        data = scalar,
        covariates = "X",
        M = M,
        alpha = alpha,
        exact = TRUE,
        n_perm = n_perm
      ),
      method = "ferman_scalar_nn"
    )
  }

  run_baba_yoshida <- function() {
    cem_data <- cem_keep(
      data = survival,
      covariates = "X",
      n_bins = n_bins,
      treat_col = "Z",
      id_col = "id"
    )
    diagnostics <- cem_match_diagnostics(cem_data, treat_col = "Z", cell_col = "cem_cell")
    by <- baba_yoshida_cem_logrank(
      data = cem_data,
      time_col = "time",
      status_col = "status",
      treat_col = "Z",
      cell_col = "cem_cell",
      id_col = "id"
    )
    by$diagnostics <- diagnostics
    by
  }

  version_a <- NULL
  if ("version_a_cem_sign" %in% methods) {
    version_a <- safe_method_result(
      cem_survival_sign_test(
        data = survival,
        covariates = "X",
        n_bins = n_bins,
        alpha = alpha,
        group_by_cell = FALSE,
        exact = TRUE,
        n_perm = n_perm
      ),
      method = "version_a_cem_sign"
    )
  }

  baba <- NULL
  if ("baba_yoshida_cem_gaussian" %in% methods) {
    if (!is.null(version_a)) {
      if (is.null(version_a$error)) {
        baba <- version_a$baba_yoshida
        baba$diagnostics <- version_a$diagnostics
      } else {
        baba <- list(method = "baba_yoshida_cem_gaussian", error = version_a$error)
      }
    } else {
      baba <- safe_method_result(
        run_baba_yoshida(),
        method = "baba_yoshida_cem_gaussian"
      )
    }
  }

  version_b <- NULL
  if ("version_b_nn_sign" %in% methods) {
    version_b <- safe_method_result(
      nn_survival_sign_test(
        data = survival,
        covariates = "X",
        M = M,
        alpha = alpha,
        exact = TRUE,
        n_perm = n_perm
      ),
      method = "version_b_nn_sign"
    )
  }

  rows <- list()
  if ("baba_yoshida_cem_gaussian" %in% methods) {
    by_p <- if (is.null(baba$error)) baba$p_gaussian_two_sided else NA_real_
    by_z <- if (is.null(baba$error)) baba$z_gaussian else NA_real_
    by_diag <- if (is.null(baba$error)) baba$diagnostics else NULL
    rows <- c(rows, list(data.frame(
      iter = iter,
      method = "baba_yoshida_cem_gaussian",
      statistic = by_z,
      p_value = by_p,
      reject = ifelse(is.na(by_p), NA, by_p <= alpha),
      n_sign_units = NA_integer_,
      n_contributions = NA_integer_,
      n_total_treated = ifelse(is.null(by_diag), n1, by_diag$n_total_treated),
      n_matched_treated = ifelse(is.null(by_diag), NA_integer_, by_diag$n_matched_treated),
      n_discarded_treated = ifelse(is.null(by_diag), NA_integer_, by_diag$n_discarded_treated),
      min_control_reuse = ifelse(is.null(by_diag) || !("min_control_reuse" %in% names(by_diag)), NA_real_, by_diag$min_control_reuse),
      mean_control_reuse = ifelse(is.null(by_diag) || !("mean_control_reuse" %in% names(by_diag)), NA_real_, by_diag$mean_control_reuse),
      max_control_reuse = ifelse(is.null(by_diag) || !("max_control_reuse" %in% names(by_diag)), NA_real_, by_diag$max_control_reuse),
      n_reused_controls = ifelse(is.null(by_diag) || !("n_reused_controls" %in% names(by_diag)), NA_integer_, by_diag$n_reused_controls),
      number_unique_controls = ifelse(is.null(by_diag) || !("number_unique_controls" %in% names(by_diag)), NA_integer_, by_diag$number_unique_controls),
      unmatched_treated = ifelse(is.null(by_diag) || !("unmatched_treated" %in% names(by_diag)), NA_integer_, by_diag$unmatched_treated),
      unmatched_controls = ifelse(is.null(by_diag) || !("unmatched_controls" %in% names(by_diag)), NA_integer_, by_diag$unmatched_controls),
      reuse_definition = ifelse(is.null(by_diag) || !("reuse_definition" %in% names(by_diag)), NA_character_, by_diag$reuse_definition),
      error = ifelse(is.null(baba$error), NA_character_, baba$error),
      stringsAsFactors = FALSE
    )))
  }

  if ("version_a_cem_sign" %in% methods) {
    rows <- c(rows, list(collect_sim_row(iter, "version_a_cem_sign", version_a, alpha)))
  }

  if ("version_b_nn_sign" %in% methods) {
    rows <- c(rows, list(collect_sim_row(iter, "version_b_nn_sign", version_b, alpha)))
  }

  if (include_scalar_ferman) {
    rows <- c(list(collect_sim_row(iter, "ferman_scalar_nn", ferman, alpha)), rows)
  }

  do.call(rbind, rows)
}

collect_sim_row <- function(iter, method, result, alpha) {
  if (!is.null(result$error)) {
    return(data.frame(
      iter = iter,
      method = method,
      statistic = NA_real_,
      p_value = NA_real_,
      reject = NA,
      n_sign_units = NA_integer_,
      n_contributions = NA_integer_,
      n_total_treated = NA_integer_,
      n_matched_treated = NA_integer_,
      n_discarded_treated = NA_integer_,
      min_control_reuse = NA_real_,
      mean_control_reuse = NA_real_,
      max_control_reuse = NA_integer_,
      n_reused_controls = NA_integer_,
      number_unique_controls = NA_integer_,
      unmatched_treated = NA_integer_,
      unmatched_controls = NA_integer_,
      reuse_definition = NA_character_,
      error = result$error,
      stringsAsFactors = FALSE
    ))
  }

  diagnostics <- result$diagnostics
  if (is.null(diagnostics)) {
    diagnostics <- data.frame(
      n_total_treated = NA_integer_,
      n_matched_treated = NA_integer_,
      n_discarded_treated = NA_integer_,
      min_control_reuse = NA_real_,
      mean_control_reuse = NA_real_,
      max_control_reuse = NA_integer_,
      n_reused_controls = NA_integer_,
      number_unique_controls = NA_integer_,
      unmatched_treated = NA_integer_,
      unmatched_controls = NA_integer_,
      reuse_definition = NA_character_
    )
  }

  data.frame(
    iter = iter,
    method = method,
    statistic = result$statistic,
    p_value = result$p_value,
    reject = result$reject,
    n_sign_units = result$n_sign_units,
    n_contributions = result$n_contributions,
    n_total_treated = diagnostics$n_total_treated,
    n_matched_treated = diagnostics$n_matched_treated,
    n_discarded_treated = diagnostics$n_discarded_treated,
    min_control_reuse = ifelse("min_control_reuse" %in% names(diagnostics), diagnostics$min_control_reuse, NA_real_),
    mean_control_reuse = ifelse("mean_control_reuse" %in% names(diagnostics), diagnostics$mean_control_reuse, NA_real_),
    max_control_reuse = ifelse("max_control_reuse" %in% names(diagnostics), diagnostics$max_control_reuse, NA_integer_),
    n_reused_controls = ifelse("n_reused_controls" %in% names(diagnostics), diagnostics$n_reused_controls, NA_integer_),
    number_unique_controls = ifelse("number_unique_controls" %in% names(diagnostics), diagnostics$number_unique_controls, NA_integer_),
    unmatched_treated = ifelse("unmatched_treated" %in% names(diagnostics), diagnostics$unmatched_treated, diagnostics$n_discarded_treated),
    unmatched_controls = ifelse("unmatched_controls" %in% names(diagnostics), diagnostics$unmatched_controls, NA_integer_),
    reuse_definition = ifelse("reuse_definition" %in% names(diagnostics), diagnostics$reuse_definition, NA_character_),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_small_simulation <- function(n_iter = 10, n1 = 5, n0 = 200, M = 2,
                                 alpha = 0.05, tau_scalar = 0,
                                 survival_hr = 1, n_bins = 4,
                                 censor_rate = 0.08, x_shift_treated = 0,
                                 seed_base = 1000,
                                 include_scalar_ferman = TRUE,
                                 n_perm = 1000) {
  rows <- lapply(seq_len(n_iter), function(iter) {
    run_one_simulation(
      iter = iter,
      n1 = n1,
      n0 = n0,
      M = M,
      alpha = alpha,
      tau_scalar = tau_scalar,
      survival_hr = survival_hr,
      n_bins = n_bins,
      censor_rate = censor_rate,
      x_shift_treated = x_shift_treated,
      seed_base = seed_base,
      include_scalar_ferman = include_scalar_ferman,
      n_perm = n_perm
    )
  })

  results <- do.call(rbind, rows)
  results$reject_numeric <- as.numeric(results$reject)
  mean_or_na <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }
  summary <- do.call(rbind, lapply(split(results, results$method), function(d) {
    data.frame(
      method = d$method[1],
      rejection_rate = mean_or_na(d$reject_numeric),
      mean_p_value = mean_or_na(d$p_value),
      mean_sign_units = mean_or_na(d$n_sign_units),
      mean_contributions = mean_or_na(d$n_contributions),
      mean_matched_treated = mean_or_na(d$n_matched_treated),
      mean_discarded_treated = mean_or_na(d$n_discarded_treated),
      mean_min_control_reuse = mean_or_na(d$min_control_reuse),
      mean_control_reuse = mean_or_na(d$mean_control_reuse),
      mean_max_control_reuse = mean_or_na(d$max_control_reuse),
      mean_reused_controls = mean_or_na(d$n_reused_controls),
      mean_number_unique_controls = mean_or_na(d$number_unique_controls),
      mean_unmatched_treated = mean_or_na(d$unmatched_treated),
      mean_unmatched_controls = mean_or_na(d$unmatched_controls),
      p_value_q10 = if (all(is.na(d$p_value))) NA_real_ else stats::quantile(d$p_value, 0.10, na.rm = TRUE),
      p_value_q50 = if (all(is.na(d$p_value))) NA_real_ else stats::quantile(d$p_value, 0.50, na.rm = TRUE),
      p_value_q90 = if (all(is.na(d$p_value))) NA_real_ else stats::quantile(d$p_value, 0.90, na.rm = TRUE),
      n_errors = sum(!is.na(d$error)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL

  list(results = results, summary = summary)
}

summarize_simulation_results <- function(results, group_cols = "method") {
  results$reject_numeric <- as.numeric(results$reject)
  mean_or_na <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }
  quantile_or_na <- function(x, prob) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    as.numeric(stats::quantile(x, prob, na.rm = TRUE))
  }

  groups <- split(results, interaction(results[, group_cols, drop = FALSE], drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(d) {
    key <- d[1, group_cols, drop = FALSE]
    data.frame(
      key,
      n_iter = length(unique(d$iter)),
      n_rows = nrow(d),
      rejection_rate = mean_or_na(d$reject_numeric),
      mean_p_value = mean_or_na(d$p_value),
      p_value_q10 = quantile_or_na(d$p_value, 0.10),
      p_value_q50 = quantile_or_na(d$p_value, 0.50),
      p_value_q90 = quantile_or_na(d$p_value, 0.90),
      mean_sign_units = mean_or_na(d$n_sign_units),
      mean_contributions = mean_or_na(d$n_contributions),
      mean_matched_treated = mean_or_na(d$n_matched_treated),
      mean_discarded_treated = mean_or_na(d$n_discarded_treated),
      mean_min_control_reuse = mean_or_na(d$min_control_reuse),
      mean_control_reuse = mean_or_na(d$mean_control_reuse),
      mean_max_control_reuse = mean_or_na(d$max_control_reuse),
      mean_reused_controls = mean_or_na(d$n_reused_controls),
      mean_number_unique_controls = mean_or_na(d$number_unique_controls),
      mean_unmatched_treated = mean_or_na(d$unmatched_treated),
      mean_unmatched_controls = mean_or_na(d$unmatched_controls),
      n_errors = sum(!is.na(d$error)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary
}

append_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    df,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    append = file.exists(path)
  )
}

cem_bins_cube_root <- function(n_total, n1 = NULL, n0 = NULL, min_bins = 4L) {
  max(min_bins, as.integer(floor(n_total^(1 / 3))))
}

survival_setting_done <- function(existing_summary, scenario, alternative_strength,
                                  survival_hr, n1, M, n_bins,
                                  x_shift_treated, censor_rate,
                                  expected_methods, n_iter) {
  if (is.null(existing_summary) || nrow(existing_summary) == 0) {
    return(FALSE)
  }
  if (!all(c("alternative_strength", "survival_hr", "n_bins") %in% names(existing_summary))) {
    return(FALSE)
  }

  rows <- existing_summary[
    existing_summary$scenario == scenario &
      existing_summary$alternative_strength == alternative_strength &
      abs(existing_summary$survival_hr - survival_hr) < 1e-12 &
      existing_summary$n1 == n1 &
      existing_summary$M == M &
      existing_summary$n_bins == n_bins &
      abs(existing_summary$x_shift_treated - x_shift_treated) < 1e-12 &
      abs(existing_summary$censor_rate - censor_rate) < 1e-12,
  ]
  if (nrow(rows) == 0) {
    return(FALSE)
  }

  all(expected_methods %in% rows$method) && all(rows$n_iter[rows$method %in% expected_methods] >= n_iter)
}

run_formal_simulation_grid <- function(n_iter = 20,
                                       n1_values = c(5, 10),
                                       n0 = 300,
                                       M_values = c(1, 2),
                                       overlap_shifts = c(0, 0.75),
                                       censor_rates = c(0.08),
                                       scenarios = c("null", "alternative"),
                                       alternative_hr_values = c(
                                         weak = 0.70,
                                         medium = 0.50,
                                         strong = 0.35
                                       ),
                                       alpha = 0.05,
                                       n_bins = 4,
                                       seed_base = 20260728,
                                       progress = FALSE,
                                       include_scalar_ferman = FALSE,
                                       n_perm = 1000,
                                       methods = c(
                                         "baba_yoshida_cem_gaussian",
                                         "version_a_cem_sign",
                                         "version_b_nn_sign"
                                       ),
                                       checkpoint_dir = NULL,
                                       checkpoint_prefix = NULL,
                                       resume = TRUE,
                                       n_cores = 1L) {
  methods <- match.arg(
    methods,
    choices = c(
      "baba_yoshida_cem_gaussian",
      "version_a_cem_sign",
      "version_b_nn_sign"
    ),
    several.ok = TRUE
  )
  all_rows <- list()
  current_summaries <- list()
  row_id <- 1
  summary_id <- 1
  grid_id <- 0
  total_grid <- length(scenarios) * length(n1_values) * length(M_values) *
    length(overlap_shifts) * length(censor_rates)
  scenario_grid <- data.frame(
    scenario = character(),
    alternative_strength = character(),
    survival_hr = numeric(),
    stringsAsFactors = FALSE
  )
  if ("null" %in% scenarios) {
    scenario_grid <- rbind(
      scenario_grid,
      data.frame(
        scenario = "null",
        alternative_strength = "null",
        survival_hr = 1,
        stringsAsFactors = FALSE
      )
    )
  }
  if ("alternative" %in% scenarios) {
    scenario_grid <- rbind(
      scenario_grid,
      data.frame(
        scenario = "alternative",
        alternative_strength = names(alternative_hr_values),
        survival_hr = as.numeric(alternative_hr_values),
        stringsAsFactors = FALSE
      )
    )
  }
  total_grid <- nrow(scenario_grid) * length(n1_values) * length(M_values) *
    length(overlap_shifts) * length(censor_rates)
  expected_methods <- c(
    if (include_scalar_ferman) "ferman_scalar_nn",
    methods
  )

  checkpoint_summary_path <- NULL
  checkpoint_raw_path <- NULL
  existing_summary <- NULL
  if (!is.null(checkpoint_dir) && !is.null(checkpoint_prefix)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    checkpoint_summary_path <- file.path(checkpoint_dir, paste0(checkpoint_prefix, "_grid_summary_all_methods.csv"))
    checkpoint_raw_path <- file.path(checkpoint_dir, paste0(checkpoint_prefix, "_raw_results.csv"))
    if (resume && file.exists(checkpoint_summary_path)) {
      existing_summary <- read.csv(checkpoint_summary_path, stringsAsFactors = FALSE)
    }
  }

  for (scenario_row in seq_len(nrow(scenario_grid))) {
    scenario <- scenario_grid$scenario[scenario_row]
    alternative_strength <- scenario_grid$alternative_strength[scenario_row]
    survival_hr <- scenario_grid$survival_hr[scenario_row]
    tau_scalar <- if (scenario == "null") 0 else 0.5

    for (n1 in n1_values) {
      for (M in M_values) {
        for (x_shift in overlap_shifts) {
          for (censor_rate in censor_rates) {
            grid_id <- grid_id + 1
            setting_n_bins <- if (is.function(n_bins)) {
              n_bins(n_total = n1 + n0, n1 = n1, n0 = n0)
            } else {
              n_bins
            }
            already_done <- survival_setting_done(
              existing_summary = existing_summary,
              scenario = scenario,
              alternative_strength = alternative_strength,
              survival_hr = survival_hr,
              n1 = n1,
              M = M,
              n_bins = setting_n_bins,
              x_shift_treated = x_shift,
              censor_rate = censor_rate,
              expected_methods = expected_methods,
              n_iter = n_iter
            )
            if (progress) {
              status <- if (already_done) "skip existing" else "run"
              cat(sprintf(
                "[survival grid %d/%d] %s: scenario=%s, strength=%s, HR=%.2f, n1=%d, M=%d, bins=%d, overlap_shift=%.2f, censor_rate=%.2f, reps=%d\n",
                grid_id, total_grid, status, scenario, alternative_strength,
                survival_hr, n1, M, setting_n_bins, x_shift, censor_rate, n_iter
              ))
              flush.console()
            }
            if (already_done) {
              next
            }

            run_iter <- function(iter) {
              one <- run_one_simulation(
                iter = iter,
                n1 = n1,
                n0 = n0,
                M = M,
                alpha = alpha,
                tau_scalar = tau_scalar,
                survival_hr = survival_hr,
                n_bins = setting_n_bins,
                censor_rate = censor_rate,
                x_shift_treated = x_shift,
                  seed_base = seed_base +
                  1000000 * scenario_row +
                  10000 * n1 +
                  1000 * M +
                  100 * round(10 * x_shift) +
                  10 * round(100 * censor_rate),
                include_scalar_ferman = include_scalar_ferman,
                n_perm = n_perm,
                methods = methods
              )
              one$scenario <- scenario
              one$alternative_strength <- alternative_strength
              one$n1 <- n1
              one$n0 <- n0
              one$M <- M
              one$tau_scalar <- tau_scalar
              one$survival_hr <- survival_hr
              one$x_shift_treated <- x_shift
              one$n_bins <- setting_n_bins
              one$censor_rate <- censor_rate
              one
            }
            if (n_cores > 1L && .Platform$OS.type != "windows") {
              setting_rows <- parallel::mclapply(
                seq_len(n_iter),
                run_iter,
                mc.cores = n_cores,
                mc.preschedule = TRUE
              )
            } else {
              setting_rows <- lapply(seq_len(n_iter), run_iter)
            }
            setting_results <- do.call(rbind, setting_rows)
            setting_summary <- summarize_simulation_results(
              setting_results,
              group_cols = c(
                "scenario", "alternative_strength", "survival_hr", "n1", "M", "n_bins",
                "x_shift_treated", "censor_rate", "method"
              )
            )
            all_rows[[row_id]] <- setting_results
            current_summaries[[summary_id]] <- setting_summary
            row_id <- row_id + 1
            summary_id <- summary_id + 1

            if (!is.null(checkpoint_summary_path)) {
              append_csv(setting_results, checkpoint_raw_path)
              append_csv(setting_summary, checkpoint_summary_path)
              existing_summary <- read.csv(checkpoint_summary_path, stringsAsFactors = FALSE)
            }
          }
        }
      }
    }
  }

  if (!is.null(checkpoint_summary_path) && file.exists(checkpoint_summary_path)) {
    summary <- read.csv(checkpoint_summary_path, stringsAsFactors = FALSE)
  } else {
            results_for_summary <- do.call(rbind, all_rows)
    summary <- summarize_simulation_results(
      results_for_summary,
      group_cols = c(
        "scenario", "alternative_strength", "survival_hr", "n1", "M", "n_bins",
        "x_shift_treated", "censor_rate", "method"
      )
    )
  }

  if (!is.null(checkpoint_raw_path) && file.exists(checkpoint_raw_path)) {
    results <- read.csv(checkpoint_raw_path, stringsAsFactors = FALSE)
  } else {
    results <- do.call(rbind, all_rows)
  }

  list(results = results, summary = summary)
}


# ============================================================
# MODULE 6: Ferman 2021 replication master
#
# This block implements the non-survival nearest-neighbor replication
# grid, including the added null/alternative tau settings.
# ============================================================

generate_ferman_dgp <- function(N1, N0, panel = "A", tau = 0, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- N1 + N0
  Z <- c(rep(1, N1), rep(0, N0))

  if (panel == "F") {
    Xmat <- matrix(stats::rnorm(n * 10), nrow = n, ncol = 10)
    colnames(Xmat) <- paste0("X", seq_len(10))
    mu1 <- 10 * sin(pi * Xmat[, 1] * Xmat[, 2]) +
      20 * (Xmat[, 3] - 0.5)^2 +
      10 * Xmat[, 4] +
      5 * Xmat[, 5] +
      Xmat[, 3] * cos(pi * Xmat[, 1] * Xmat[, 2])
    tau_x <- Xmat[, 3] * cos(pi * Xmat[, 1] * Xmat[, 2])
    mu0 <- mu1 - tau_x
    eps1 <- stats::rnorm(n)
  } else {
    X <- stats::rnorm(n)
    Xmat <- matrix(X, ncol = 1)
    colnames(Xmat) <- "X"
    mu0 <- rep(0, n)

    if (panel == "A") {
      mu1 <- X
      eps1 <- stats::rnorm(n)
    } else if (panel == "B") {
      mu1 <- (stats::qchisq(stats::pnorm(X), df = 8) - 8) / sqrt(16)
      eps1 <- stats::rnorm(n)
    } else if (panel == "C") {
      mu1 <- (stats::qchisq(stats::pnorm(X), df = 1) - 1) / sqrt(2)
      eps1 <- stats::rnorm(n)
    } else if (panel == "D") {
      mu1 <- (stats::qchisq(stats::pnorm(X), df = 1) - 1) / sqrt(2)
      eps1 <- (stats::rchisq(n, df = 1) - 1) / sqrt(2)
    } else if (panel == "E") {
      mu1 <- (stats::qchisq(stats::pnorm(X), df = 1) - 1) / sqrt(2)
      eps1 <- 2 * (stats::rchisq(n, df = 1) - 1) / sqrt(2)
    } else {
      stop("panel must be one of A, B, C, D, E, F.")
    }
  }

  eps0 <- stats::rnorm(n)
  Y0 <- mu0 + eps0
  Y1 <- mu1 + eps1
  Y <- ifelse(Z == 1, Y1 + tau, Y0)

  data.frame(
    id = paste0(ifelse(Z == 1, "t", "c"), seq_len(n)),
    Z = Z,
    Xmat,
    Y0 = Y0,
    Y1 = Y1,
    Y = Y,
    noise_0 = eps0,
    noise_1 = eps1,
    check.names = FALSE
  )
}

run_one_ferman_replication <- function(iter = 1, N1 = 10, N0 = 1000,
                                       M = 4, panel = "A", tau = 0, tau0 = 0,
                                       scenario = ifelse(tau == 0, "null", "alternative"),
                                       alternative_strength = ifelse(tau == 0, "null", "standard"),
                                       alpha = 0.05, n_perm = 1000,
                                       seed_base = 7000) {
  dat <- generate_ferman_dgp(
    N1 = N1,
    N0 = N0,
    panel = panel,
    tau = tau,
    seed = seed_base +
      100000 * match(panel, c("A", "B", "C", "D", "E", "F")) +
      1000 * N1 +
      100 * M +
      10000 * as.integer(round(abs(tau) * 100)) +
      iter
  )
  covariates <- grep("^X", names(dat), value = TRUE)
  res <- ferman_scalar_nn_test(
    data = dat,
    covariates = covariates,
    outcome_col = "Y",
    M = M,
    tau0 = tau0,
    alpha = alpha,
    exact = TRUE,
    n_perm = n_perm,
    seed = seed_base + iter
  )

  data.frame(
    iter = iter,
    scenario = scenario,
    alternative_strength = alternative_strength,
    tau = tau,
    panel = panel,
    N1 = N1,
    N0 = N0,
    M = M,
    statistic = res$statistic,
    p_value = res$p_value,
    reject = res$reject,
    n_sign_units = res$n_sign_units,
    n_contributions = res$n_contributions,
    n_matched_treated = res$diagnostics$n_matched_treated,
    max_control_reuse = res$diagnostics$max_control_reuse,
    n_reused_controls = res$diagnostics$n_reused_controls,
    stringsAsFactors = FALSE
  )
}

summarize_ferman_replication <- function(results) {
  results$reject_numeric <- as.numeric(results$reject)
  mean_or_na <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }

  if (!("scenario" %in% names(results))) {
    results$scenario <- "null"
  }
  if (!("alternative_strength" %in% names(results))) {
    results$alternative_strength <- "null"
  }
  if (!("tau" %in% names(results))) {
    results$tau <- 0
  }

  group_cols <- c("scenario", "alternative_strength", "tau", "panel", "N1", "M")
  groups <- split(results, interaction(results[, group_cols, drop = FALSE], drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(d) {
    data.frame(
      scenario = d$scenario[1],
      alternative_strength = d$alternative_strength[1],
      tau = d$tau[1],
      panel = d$panel[1],
      N1 = d$N1[1],
      M = d$M[1],
      n_iter = length(unique(d$iter)),
      rejection_rate = mean_or_na(d$reject_numeric),
      mean_p_value = mean_or_na(d$p_value),
      mean_sign_units = mean_or_na(d$n_sign_units),
      mean_max_control_reuse = mean_or_na(d$max_control_reuse),
      mean_reused_controls = mean_or_na(d$n_reused_controls),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary[order(summary$scenario, summary$alternative_strength, summary$panel, summary$N1, summary$M), ]
}

ferman_setting_done <- function(existing_summary, panel, N1, M, n_iter,
                                scenario = "null", alternative_strength = "null",
                                tau = 0) {
  if (is.null(existing_summary) || nrow(existing_summary) == 0) {
    return(FALSE)
  }
  if (!("scenario" %in% names(existing_summary))) {
    existing_summary$scenario <- "null"
  }
  if (!("alternative_strength" %in% names(existing_summary))) {
    existing_summary$alternative_strength <- "null"
  }
  if (!("tau" %in% names(existing_summary))) {
    existing_summary$tau <- 0
  }

  rows <- existing_summary[
    existing_summary$scenario == scenario &
      existing_summary$alternative_strength == alternative_strength &
      abs(existing_summary$tau - tau) < 1e-12 &
      existing_summary$panel == panel &
      existing_summary$N1 == N1 &
      existing_summary$M == M,
  ]
  nrow(rows) > 0 && all(rows$n_iter >= n_iter)
}

run_ferman_replication_grid <- function(n_iter = 20,
                                        panels = c("A", "B", "C", "D", "E"),
                                        N1_values = c(5, 10, 25, 50),
                                        M_values = c(1, 4, 10),
                                        N0 = 1000,
                                        tau_values = c(null = 0),
                                        tau0 = 0,
                                        alpha = 0.05,
                                        n_perm = 1000,
                                        seed_base = 7000,
                                        progress = FALSE,
                                        checkpoint_dir = NULL,
                                        checkpoint_prefix = NULL,
                                        resume = TRUE) {
  rows <- list()
  current_summaries <- list()
  row_id <- 1
  summary_id <- 1
  grid_id <- 0
  if (is.null(names(tau_values)) || any(names(tau_values) == "")) {
    names(tau_values) <- ifelse(tau_values == 0, "null", paste0("tau_", tau_values))
  }
  total_grid <- length(tau_values) * length(panels) * length(N1_values) * length(M_values)
  checkpoint_summary_path <- NULL
  checkpoint_raw_path <- NULL
  existing_summary <- NULL
  if (!is.null(checkpoint_dir) && !is.null(checkpoint_prefix)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    checkpoint_summary_path <- file.path(checkpoint_dir, paste0(checkpoint_prefix, "_grid_summary.csv"))
    checkpoint_raw_path <- file.path(checkpoint_dir, paste0(checkpoint_prefix, "_raw_results.csv"))
    if (resume && file.exists(checkpoint_summary_path)) {
      existing_summary <- read.csv(checkpoint_summary_path, stringsAsFactors = FALSE)
    }
  }

  for (tau_name in names(tau_values)) {
    tau <- unname(tau_values[[tau_name]])
    scenario <- if (tau == 0) "null" else "alternative"
    alternative_strength <- if (tau == 0) "null" else tau_name
    for (panel in panels) {
      for (N1 in N1_values) {
        for (M in M_values) {
          grid_id <- grid_id + 1
          already_done <- ferman_setting_done(
            existing_summary,
            panel = panel,
            N1 = N1,
            M = M,
            n_iter = n_iter,
            scenario = scenario,
            alternative_strength = alternative_strength,
            tau = tau
          )
        if (progress) {
          status <- if (already_done) "skip existing" else "run"
          cat(sprintf(
            "[ferman grid %d/%d] %s: scenario=%s, tau=%.3f, panel=%s, N1=%d, M=%d, reps=%d\n",
            grid_id, total_grid, status, scenario, tau, panel, N1, M, n_iter
          ))
          flush.console()
        }
        if (already_done) {
          next
        }
        setting_rows <- vector("list", n_iter)
        for (iter in seq_len(n_iter)) {
          setting_rows[[iter]] <- run_one_ferman_replication(
            iter = iter,
            N1 = N1,
            N0 = N0,
            M = M,
            panel = panel,
            tau = tau,
            scenario = scenario,
            alternative_strength = alternative_strength,
            tau0 = tau0,
            alpha = alpha,
            n_perm = n_perm,
            seed_base = seed_base
          )
        }
        setting_results <- do.call(rbind, setting_rows)
        setting_summary <- summarize_ferman_replication(setting_results)
        rows[[row_id]] <- setting_results
        current_summaries[[summary_id]] <- setting_summary
        row_id <- row_id + 1
        summary_id <- summary_id + 1

        if (!is.null(checkpoint_summary_path)) {
          append_csv(setting_results, checkpoint_raw_path)
          append_csv(setting_summary, checkpoint_summary_path)
          existing_summary <- read.csv(checkpoint_summary_path, stringsAsFactors = FALSE)
        }
      }
    }
    }
  }

  if (!is.null(checkpoint_raw_path) && file.exists(checkpoint_raw_path)) {
    results <- read.csv(checkpoint_raw_path, stringsAsFactors = FALSE)
  } else {
    results <- do.call(rbind, rows)
  }
  if (!is.null(checkpoint_summary_path) && file.exists(checkpoint_summary_path)) {
    summary <- read.csv(checkpoint_summary_path, stringsAsFactors = FALSE)
  } else {
    summary <- summarize_ferman_replication(results)
  }

  list(results = results, summary = summary)
}

compare_to_old_ferman_outputs <- function(our_summary,
                                          old_path = "results/old_attempts_ferman_saved_results_summary.csv",
                                          old_source = "results_table.rds") {
  old <- read.csv(old_path, stringsAsFactors = FALSE)
  old <- old[old$source_file == old_source, ]
  names(old)[names(old) == "rejection_rate"] <- "old_attempt_rejection_rate"
  merged <- merge(
    our_summary,
    old[, c("panel", "N1", "M", "old_attempt_rejection_rate") ],
    by = c("panel", "N1", "M"),
    all.x = TRUE
  )
  merged$rejection_rate_diff <- merged$rejection_rate - merged$old_attempt_rejection_rate
  merged[order(merged$panel, merged$N1, merged$M), ]
}
