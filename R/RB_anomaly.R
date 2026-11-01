# Shared record engine: gauges x chronological years, no detrending.
#
# Reading order / call flow:
#   validation -> rb_scan_records() -> rb_aggregate_records() -> rb_anomaly_ci()
# Public vector and regional functions reuse this path. Legacy CI wrappers
# adapt the output fields while retaining their original statistical targets.
# Internal helpers use ordinary comments rather than exported help topics.

# ---- Input validation and preparation ------------------------------------

# Validate one gauge-by-year dataset and return a numeric matrix.
# Accept NA/NaN as missing observations, but reject infinite values and text.
# Data-frame columns must already be numeric: conversion must not silently
# turn factor levels or character labels into hydrological observations.
rb_validate_matrix <- function(x) {
  # Reject unsupported containers and zero-sized dimensions before conversion.
  if (!(is.matrix(x) || is.data.frame(x)) || any(dim(x) == 0L)) {
    stop("Input must be a non-empty numeric matrix or data frame.")
  }
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop("All input columns must be numeric.")
  }
  # Preserve matrix dimensions and row/column names for detailed outputs.
  x <- as.matrix(x)
  if (!is.numeric(x) || any(is.infinite(x))) {
    stop("Input must contain numeric, finite observations or NA.")
  }
  x
}

# Validate the calendar associated with n input columns; return it invisibly.
# Never reorder data here: column order defines each gauge's record history.
# regular = TRUE adds the equal-spacing requirement of the SSA smoother.
rb_validate_year <- function(year, n, regular = FALSE) {
  # A unique increasing calendar prevents accidental reversal or duplicate years.
  if (!is.numeric(year) || length(year) != n || any(!is.finite(year)) ||
      any(diff(year) <= 0)) {
    stop("Length of year must match columns; years must be finite and strictly increasing.")
  }
  # Compare gaps with numerical tolerance rather than requiring exact equality.
  if (regular && length(year) > 2L &&
      !isTRUE(all.equal(diff(year), rep(diff(year)[1], length(year) - 1L)))) {
    stop("SSA requires regularly spaced years; supply missing years explicitly.")
  }
  invisible(year)
}

# Validate a non-empty list of aligned datasets and their shared calendar.
# Return the same list structure with numeric matrices ready for computation.
# Shape checks do not establish matching gauge identities: callers must align
# rows and columns across months/variables before supplying the list.
rb_prepare_matrices <- function(list_df, year) {
  # A data frame is technically an R list, but is one dataset, not a list of them.
  if (!is.list(list_df) || is.data.frame(list_df) || !length(list_df)) {
    stop("list_df must be a non-empty list of matrices or data frames.")
  }
  # Convert once outside Monte Carlo loops; preserve list names.
  out <- lapply(list_df, rb_validate_matrix)
  # All list elements must share the first element's gauge and year dimensions.
  if (!all(vapply(out, function(x) identical(dim(x), dim(out[[1]])), logical(1)))) {
    stop("Elements in list_df have different shapes.")
  }
  # One calendar check suffices once all matrix dimensions agree.
  rb_validate_year(year, ncol(out[[1]]))
  out
}

# Guard a scalar loop/allocation parameter such as runs.
# label supplies the parameter name in errors; successful calls only validate.
# Integer-valued doubles (e.g. 1000) are allowed within R's integer range.
rb_validate_positive_integer <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(label, " must be a positive integer.")
  }
}

# ---- Record engine and count aggregation ---------------------------------

# Scan a prevalidated numeric gauge-by-year matrix in the supplied order.
# perm must contain every column index exactly once; internal callers generate
# it with seq_len() or sample.int(), so validation is not repeated here.
# detailed = TRUE returns gauge-level Robs/R1n matrices; FALSE returns only
# year-position totals Robs/R1n/n_eligible, reducing simulation memory use.
# Output positions follow perm, not original calendar labels. The detailed
# mode is used publicly only with original ordering (and original dimnames).
# O(gauges * years). Index input directly instead of copying shuffled matrices.
# Compare current values before updating running maxima and valid counts.
rb_scan_records <- function(x, perm = seq_len(ncol(x)), detailed = FALSE) {
  p <- nrow(x)
  n <- ncol(x)
  # Each gauge keeps its own history. -Inf initializes the maximum; past_n = 0
  # prevents the first valid observation from being counted as a record.
  past_max <- rep(-Inf, p)
  past_n <- integer(p)
  # Zero is the sum over an empty eligible set; its anomaly will become NA.
  observed <- expected <- numeric(n)
  eligible <- integer(n)
  # Detailed entries remain NA unless a valid comparison with history exists.
  if (detailed) Robs <- R1n <- matrix(NA_real_, p, n, dimnames = dimnames(x))
  for (j in seq_len(n)) {
    # j is the new sequence position; perm[j] is the original input column.
    now <- x[, perm[j]]
    valid <- !is.na(now)
    # The same eligibility mask must govern both observed and expected counts.
    use <- valid & past_n > 0L
    # Strict comparison treats ties as non-records; expected probabilities use
    # valid history INCLUDING the current value, hence past_n + 1.
    obs <- as.numeric(now[use] > past_max[use])
    exp <- 1 / (past_n[use] + 1)
    # Sum counts, not percentages, across gauges at this sequence position.
    observed[j] <- sum(obs)
    expected[j] <- sum(exp)
    eligible[j] <- sum(use)
    if (detailed) {
      # Populate only eligible gauge cells; retain NA for missing/first values.
      Robs[use, j] <- obs
      R1n[use, j] <- exp
    }
    # Update AFTER comparison, including each gauge's first valid observation.
    # Missing observations change neither the maximum nor the history count.
    past_max[valid] <- pmax(past_max[valid], now[valid])
    past_n[valid] <- past_n[valid] + 1L
  }
  # Both output modes derive from exactly the same eligibility and comparisons.
  if (detailed) return(list(Robs = Robs, R1n = R1n))
  list(Robs = observed, R1n = expected, n_eligible = eligible)
}

# Convert equal-length non-missing count vectors to percentage anomalies.
# Internal callers provide nonnegative expectations. A zero expectation means
# there was no eligible comparison, so leave that position NA, not zero percent.
rb_compute_anomaly <- function(observed, expected) {
  out <- rep(NA_real_, length(expected))
  use <- expected > 0
  # Retain full precision; rounding belongs in presentation, not simulation.
  out[use] <- 100 * (observed[use] - expected[use]) / expected[use]
  out
}

# Pool prevalidated, aligned matrices (e.g. months within one region).
# Return one row per sequence position: observed/expected counts, percent
# anomaly, and number of eligible gauge-series contributions. perm is shared
# across every element; each element's gauges still retain separate histories.
rb_aggregate_records <- function(matrices, perm = seq_len(ncol(matrices[[1]]))) {
  n <- ncol(matrices[[1]])
  observed <- expected <- numeric(n)
  eligible <- integer(n)
  for (x in matrices) {
    # Start a fresh record history for each dataset and sum only its outputs.
    out <- rb_scan_records(x, perm)
    observed <- observed + out$Robs
    expected <- expected + out$R1n
    eligible <- eligible + out$n_eligible
  }
  # Normalize pooled counts once; averaging individual anomalies is different.
  data.frame(Robs = observed, R1n = expected,
             Ranom = rb_compute_anomaly(observed, expected), n_eligible = eligible)
}

# ---- Public record and regional interfaces -------------------------------

#' Calculate Station-Level Record Statistics
#'
#' Strict records exceed all previous valid observations. Ties are not records.
#' Missing values do not advance history. Missing and first valid observations
#' return NA in both matrices. The continuous iid expectation is 1/n, where n
#' includes the current valid observation. Ties can invalidate this analytical
#' expectation. No detrending is applied.
#' @param df Numeric matrix or data frame: gauges in rows, ordered years in columns.
#' @return List with Robs (0/1 indicators) and R1n (probabilities), retaining
#'   input dimensions and dimnames.
#' @export
#' @examples
#' rb_record_matrix(rbind(c(1, 2, 2, 4), c(2, NA, 3, 1)))
# Validate once, then request gauge-level matrices rather than regional totals.
rb_record_matrix <- function(df) rb_scan_records(rb_validate_matrix(df), detailed = TRUE)

# Adapt a single chronological series to the matrix engine.
# Return the one-row Robs/R1n matrix list; the public wrappers flatten it below.
# Preserve the older vector API's requirement of at least two valid values,
# although the matrix API allows gauges with insufficient history.
rb_record_vector <- function(x) {
  if (!is.numeric(x) || !is.null(dim(x)) || any(is.infinite(x))) {
    stop("x must be a numeric vector containing finite observations or NA.")
  }
  if (sum(!is.na(x)) < 2L) stop("At least two valid observations are required.")
  # A series is one gauge (one row), not a column of separate gauges.
  rb_record_matrix(matrix(x, nrow = 1))
}

#' Identify Record-Breaking Observations
#' @param x Numeric vector with at least two valid observations, in time order.
#' @return Numeric vector: 1 for strict records, 0 for non-records, NA for
#'   missing observations and the first valid observation.
#' @export
# Extract indicators as a plain vector, preserving sequence length and NA cells.
rb_record_indicator <- function(x) as.numeric(rb_record_vector(x)$Robs)

#' Expected Record Probabilities
#' @inheritParams rb_record_indicator
#' @return Numeric vector of 1/n probabilities; missing observations and the
#'   first valid observation return NA. Missing values do not advance n.
#' @export
# Extract expected probabilities from the same engine and eligibility rule.
rb_record_probability <- function(x) as.numeric(rb_record_vector(x)$R1n)

#' Aggregate Record Counts Within a Region
#' @inheritParams rb_record_matrix
#' @param year Numeric, finite, strictly increasing years, one per column.
#' @param plot Logical; print the anomaly plot when TRUE.
#' @return Data frame with year, R1n, Robs, Ranom (percent), and n_eligible.
#'   Ineligible years have zero counts and NA anomaly. Calculations retain
#'   full precision. An eligible observation is valid and has prior history.
#' @export
#' @examples
#' rb_anomaly_region(rbind(1:4, 4:1), 2001:2004)
rb_anomaly_region <- function(df, year, plot = FALSE) {
  # Wrap one dataset so single- and multi-dataset paths share validation.
  matrices <- rb_prepare_matrices(list(df), year)
  out <- rb_aggregate_records(matrices)
  # Keep the historical year/R1n/Robs/Ranom column order; append eligibility.
  out <- data.frame(year = year, out[c("R1n", "Robs", "Ranom", "n_eligible")])
  # Optional plotting is a side effect; the numerical result is always returned.
  if (plot) rb_plot_anomaly(out$Ranom, year)
  out
}

#' Aggregate Counts Across Aligned Datasets
#' @param list_df Non-empty list of numeric gauge-by-year matrices or data
#'   frames of identical dimensions. Elements represent aligned months or
#'   comparable series within ONE regional window, not independent regions.
#' @inheritParams rb_anomaly_region
#' @return List of per-dataset regional results and out_Ranom, containing
#'   time, Robs, R1n, Ranom (percent), and n_eligible. Counts are pooled before
#'   normalization; n_eligible counts gauge-series contributions.
#' @export
rb_anomaly_region_list <- function(list_df, year, plot = FALSE) {
  matrices <- rb_prepare_matrices(list_df, year)
  # Retain each dataset's result for inspection as well as the pooled result.
  out <- lapply(matrices, function(x) rb_anomaly_region(x, year))
  # Prevent the pooled result from overwriting a user-named dataset result.
  if ("out_Ranom" %in% names(out)) stop("out_Ranom is a reserved list name.")
  # Recompute pooled counts rather than averaging per-dataset percentages.
  out$out_Ranom <- data.frame(time = year, rb_aggregate_records(matrices))
  if (plot) rb_plot_anomaly(out$out_Ranom$Ranom, year)
  out
}

# ---- Historical internal adapters (not exported) -------------------------

# Convert a non-empty list of shaped datasets without a supplied calendar.
# A synthetic 1..n calendar validates dimensions only, not real date alignment.
# Return the numeric matrix list; the actual year check belongs to the caller.
rb_prepare_matrix_list <- function(list_df) {
  rb_prepare_matrices(list_df, seq_len(ncol(list_df[[1]])))
}

# Return Robs/R1n/n_eligible totals for a prepared numeric matrix.
# NULL means original order; a supplied perm is trusted as a full permutation.
rb_record_sums <- function(x, perm = NULL) {
  if (is.null(perm)) perm <- seq_len(ncol(x))
  rb_scan_records(x, perm)
}

# Return pooled counts and anomaly for an already prepared matrix list.
# Despite its historical name, this helper returns count/eligibility columns
# too. It intentionally omits dates, smoothing, and plot construction.
rb_anomaly_from_matrices <- function(list_mat, perm = NULL) {
  if (is.null(perm)) perm <- seq_len(ncol(list_mat[[1]]))
  rb_aggregate_records(list_mat, perm)
}

# ---- Statistic selection and null-ensemble summaries ---------------------

# Transform a raw anomaly vector into the requested comparison statistic.
# Return a vector on the same time axis. The public caller validates statistic;
# each smoother validates its own parameters and missing-value requirements.
# Calling this for both observed and simulated data keeps the bands comparable.
rb_apply_statistic <- function(x, year, statistic, L, groups, span) {
  # Raw mode bypasses smoothing. SSA uses sequence spacing; LOESS uses year.
  switch(statistic,
    raw = x,
    ssa = rb_ssa_trend(x, L = L, groups = groups),
    loess = rb_loess_trend(data.frame(time = year, value = x), span)
  )
}

# Summarize a simulations-by-time matrix with pointwise central null bands.
# Return a time-row data frame with lower/upper 90% and 95% bounds. Quantiles
# use available runs at each position; no available runs means NA bounds.
rb_null_quantiles <- function(simulations) {
  # Work down columns: runs are replicates, not consecutive time observations.
  q <- vapply(seq_len(ncol(simulations)), function(j) {
    x <- simulations[, j]
    if (all(is.na(x))) return(rep(NA_real_, 4))
    # R's default type-7 quantiles; suppress percentile names for stable fields.
    stats::quantile(x, c(0.05, 0.95, 0.025, 0.975), na.rm = TRUE, names = FALSE)
  }, numeric(4))
  # vapply produces four bounds x time; transpose to one output row per time.
  out <- as.data.frame(t(q))
  names(out) <- c("ci_lower_90", "ci_upper_90", "ci_lower_95", "ci_upper_95")
  out
}

#' Simulate Pointwise Record-Anomaly Null Reference Bands
#'
#' Permute complete year columns synchronously across all gauges and list
#' elements within one regional window. Missing entries move with their columns;
#' their original dates are not frozen. Recompute valid history each run.
#' Apply the SAME statistic to observations and each simulation before taking
#' quantiles (R quantile type 7). Bands are pointwise null-reference intervals,
#' not uncertainty intervals around the observed trend or simultaneous bands.
#'
#' The null treats years as exchangeable; temporal persistence is destroyed.
#' Continuous iid observations give the 1/n analytical baseline. Ties,
#' informative missingness and temporal dependence require separate assessment.
#' This function does not select regional windows, test iid, or shuffle blocks.
#' SSA uses regular time spacing, preserves missing edges, and rejects interior
#' gaps. L = 15 follows the paper's window; groups = list(trend = 1) is a package
#' choice, not a verified author component-selection recipe. LOESS is optional.
#' @inheritParams rb_anomaly_region_list
#' @param statistic Target of the bands: SSA trend (default), raw anomaly, or
#'   LOESS trend. Observations and simulations use the same target.
#' @param runs Positive integer simulation count; default 1000. The paper used
#'   10000; assess Monte Carlo precision for substantive analyses.
#' @param seed Nonnegative integer seed, or NULL to use the caller's RNG stream.
#'   A supplied seed is local: the caller's RNG state is restored afterward.
#' @param L SSA window length, at least 2 and smaller than usable series length.
#' @param groups Named SSA eigentriple groups containing trend.
#' @param span Positive LOESS smoothing span.
#' @return List with data (time, Robs, R1n, Ranom, n_eligible, statistic_value,
#'   and ci_lower/upper_90/95), and settings recording the statistic, simulation
#'   count, seed, smoothing parameters, permutation and missingness policies,
#'   interval probabilities, and quantile type. Values and bands are percentages.
#' @export
#' @examples
#' x <- rbind(1:20, 20:1)
#' fit <- rb_anomaly_ci(list(x), 2001:2020, statistic = "raw", runs = 10)
#' fit$data
rb_anomaly_ci <- function(list_df, year, statistic = c("ssa", "raw", "loess"),
                          runs = 1000, seed = 1, L = 15,
                          groups = list(trend = 1), span = 0.8) {
  # Resolve the selected statistic and prepare inputs before allocating runs.
  statistic <- match.arg(statistic)
  matrices <- rb_prepare_matrices(list_df, year)
  rb_validate_positive_integer(runs, "runs")
  # Only SSA assumes equally spaced input time steps; raw counting uses order.
  rb_validate_year(year, ncol(matrices[[1]]), regular = statistic == "ssa")
  if (!is.null(seed)) {
    # NULL deliberately uses and advances the caller's existing RNG stream.
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      stop("seed must be a nonnegative integer or NULL.")
    }
    # Save RNG state BEFORE any fitting or shuffling, including its absence.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    # Restore on success AND error; do not leave a new RNG state behind when
    # the caller originally had none. add = TRUE preserves any other cleanup.
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    # Establish a repeatable local stream for the entire analysis.
    set.seed(seed)
  }
  # Compute the observed anomaly on original, undetrended observations first.
  # Apply the target smoother now so invalid observed-series settings fail early.
  observed <- rb_aggregate_records(matrices)
  observed_statistic <- rb_apply_statistic(observed$Ranom, year, statistic, L, groups, span)
  n <- length(year)
  # Store only the selected statistic per run/time, not gauge-level records.
  simulated <- matrix(NA_real_, runs, n)
  for (r in seq_len(runs)) {
    # One permutation for every gauge and aligned dataset preserves their joint
    # same-year configuration while destroying chronological dependence.
    perm <- sample.int(n)
    # Recompute record histories and 1/n expectations in the shuffled order;
    # original-calendar expectations cannot be reused after moving missing data.
    values <- rb_aggregate_records(matrices, perm)$Ranom
    # Smooth each realization BEFORE computing ensemble quantiles. Report the
    # failing run rather than silently dropping an unsuitable SSA/LOESS fit.
    simulated[r, ] <- tryCatch(
      rb_apply_statistic(values, year, statistic, L, groups, span),
      error = function(e) stop("Simulation ", r, ": ", conditionMessage(e), call. = FALSE)
    )
  }
  # Attach actual year labels to sequence positions only after simulation.
  # Keep raw anomalies beside the comparison statistic, plus explicit settings
  # so downstream plots/reports can identify what the null bands describe.
  list(
    data = data.frame(time = year, observed, statistic_value = observed_statistic,
                      rb_null_quantiles(simulated)),
    settings = list(statistic = statistic, runs = runs, seed = seed,
      L = if (statistic == "ssa") L else NULL,
      groups = if (statistic == "ssa") groups else NULL,
      span = if (statistic == "loess") span else NULL,
      permutation = "shared year columns within one regional window",
      missingness = "moves with columns; valid history recomputed",
      intervals = "pointwise null-reference bands",
      probabilities = c(0.05, 0.95, 0.025, 0.975), quantile_type = 7L)
  )
}

# ---- Compatibility interfaces -------------------------------------------

# Adapt a new result (list with data) to the historical CI data-frame columns.
# loess = TRUE adds an observed descriptive LOESS line to RAW anomaly bands;
# it does not alter those bands. Otherwise copy the selected SSA/LOESS statistic
# into the legacy loess_trend field, whose name is retained for compatibility.
rb_as_legacy_ci <- function(fit, loess = FALSE) {
  x <- fit$data
  trend <- if (loess) rb_loess_trend(data.frame(time = x$time, value = x$Ranom))
    else x$statistic_value
  # Preserve original column names/order; omit new count and settings fields.
  data.frame(time = x$time, value = x$Ranom, loess_trend = trend,
    x[c("ci_lower_90", "ci_upper_90", "ci_lower_95", "ci_upper_95")])
}

#' Legacy Raw-Anomaly Bands for Aligned Datasets
#' @inheritParams rb_anomaly_ci
#' @return Data frame with time, value, loess_trend, and 90/95 percent bands.
#'   Bands apply to raw anomalies, NOT the accompanying observed LOESS trend.
#' @keywords internal
rb_legacy_raw_ci_list <- function(list_df, year, seed = 1, runs = 1000) {
  # Historical target: raw bands, with a separate observed LOESS summary line.
  rb_as_legacy_ci(rb_anomaly_ci(list_df, year, "raw", runs, seed), loess = TRUE)
}

#' Legacy Raw-Anomaly Bands for One Dataset
#' @inheritParams rb_record_matrix
#' @inheritParams rb_anomaly_ci
#' @return Data frame as in [Func_RB_region_list_CI()]. Uses the caller's RNG.
#' @keywords internal
rb_legacy_raw_ci <- function(df, year, runs = 1000) {
  # Adapt one matrix and preserve caller-controlled set.seed() behavior.
  rb_legacy_raw_ci_list(list(df), year, seed = NULL, runs = runs)
}

#' Legacy SSA-Trend Null Bands
#' @inheritParams rb_anomaly_ci
#' @return Data frame with time, value, loess_trend, and 90/95 percent bands.
#'   For compatibility loess_trend contains SSA, and bands apply to that SSA
#'   statistic. Prefer [rb_anomaly_ci()] for explicit target and settings.
#' @keywords internal
rb_legacy_ssa_ci <- function(list_df, year, seed = 1, runs = 1000,
                                        L = 15, groups = list(trend = 1)) {
  # "Fast" historically selected SSA; it is not interchangeable with raw bands.
  # All modern wrappers now benefit from the same aggregate-only engine.
  rb_as_legacy_ci(rb_anomaly_ci(list_df, year, "ssa", runs, seed, L, groups))
}

#' Legacy LOESS-Trend Null Bands
#' @inheritParams rb_record_matrix
#' @inheritParams rb_anomaly_ci
#' @return Data frame with time, value, loess_trend, ci_lower and ci_upper
#'   (90 percent LOESS-trend null bands). Uses the caller's RNG.
#' @keywords internal
rb_legacy_loess_ci <- function(df, year, runs = 1000) {
  # Request matched LOESS bands, then retain only the legacy 90% output fields.
  out <- rb_as_legacy_ci(rb_anomaly_ci(list(df), year, "loess", runs, seed = NULL))
  # The adapter guarantees lower/upper 90% bounds occupy columns four and five.
  names(out)[4:5] <- c("ci_lower", "ci_upper")
  out[1:5]
}

#' Deprecated Calendar-Reindexed Anomaly Bands
#'
#' Retained only for legacy analyses. Reassigns shuffled sequence-position
#' anomalies to the original observation years. This is NOT the sequence-position
#' null reference for the documented method. Use [RB_Ranom_CI()] or
#' [rb_anomaly_ci()] instead; their results intentionally differ.
#' @inheritParams rb_record_matrix
#' @inheritParams rb_anomaly_ci
#' @return Legacy data frame with time, value, loess_trend and 90/95 percent bands.
#' @keywords internal
rb_legacy_reindexed_ci <- function(df, year, runs = 1000) {
  # Warn but retain the historical algorithm; silently substituting standard
  # null bands would change the meaning of an existing user's calculation.
  .Deprecated("rb_anomaly_ci", msg = paste(
    "RB_Ranom_CI_index is deprecated: calendar-reindexed bands are not",
    "sequence-position null bands. Use rb_anomaly_ci(statistic = 'raw')."))
  matrices <- rb_prepare_matrices(list(df), year)
  rb_validate_positive_integer(runs, "runs")
  n <- length(year)
  simulated <- matrix(NA_real_, runs, n)
  for (r in seq_len(runs)) {
    perm <- sample.int(n)
    # LEGACY DIFFERENCE: place sequence-position j's anomaly into original
    # column perm[j]. Standard simulation stores it in column j instead.
    # Mixing record-history positions this way is why this method is deprecated.
    simulated[r, perm] <- rb_aggregate_records(matrices, perm)$Ranom
  }
  # Observed values stay in chronological order; the adapter adds their LOESS
  # summary and retains the historical calendar-reindexed quantile columns.
  observed <- rb_aggregate_records(matrices)
  rb_as_legacy_ci(list(data = data.frame(time = year, observed,
                                      rb_null_quantiles(simulated))), loess = TRUE)
}

# ---- Compatibility aliases ----------------------------------------------
# Older scripts can keep their names and named arguments. These aliases
# refer to the same implementation; use the rb_* API for new analyses.

#' @rdname rb_record_matrix
#' @details `Func_RB_matrix()` is a compatibility alias. Use `rb_record_matrix()`
#'   in new code; arguments and returned values are unchanged.
#' @export
Func_RB_matrix <- rb_record_matrix

#' @rdname rb_record_indicator
#' @details `Func_Ranom()` is a compatibility alias. Use `rb_record_indicator()`
#'   in new code; arguments and returned values are unchanged.
#' @export
Func_Ranom <- rb_record_indicator

#' @rdname rb_record_probability
#' @details `Func_Ranom_1N()` is a compatibility alias. Use `rb_record_probability()`
#'   in new code; arguments and returned values are unchanged.
#' @export
Func_Ranom_1N <- rb_record_probability

#' @rdname rb_anomaly_region
#' @details `Func_RB_region()` is a compatibility alias. Use `rb_anomaly_region()`
#'   in new code; arguments and returned values are unchanged.
#' @export
Func_RB_region <- rb_anomaly_region

#' @rdname rb_anomaly_region_list
#' @details `Func_RB_region_list()` is a compatibility alias. Use `rb_anomaly_region_list()`
#'   in new code; arguments and returned values are unchanged.
#' @export
Func_RB_region_list <- rb_anomaly_region_list

#' @rdname rb_legacy_raw_ci_list
#' @export
Func_RB_region_list_CI <- rb_legacy_raw_ci_list

#' @rdname rb_legacy_raw_ci
#' @export
RB_Ranom_CI <- rb_legacy_raw_ci

#' @rdname rb_legacy_ssa_ci
#' @export
Func_RB_region_list_CI_fast <- rb_legacy_ssa_ci

#' @rdname rb_legacy_loess_ci
#' @export
LOESS_Trend_CI <- rb_legacy_loess_ci

#' @rdname rb_legacy_reindexed_ci
#' @export
RB_Ranom_CI_index <- rb_legacy_reindexed_ci
