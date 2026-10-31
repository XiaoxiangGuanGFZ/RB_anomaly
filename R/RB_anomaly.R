# Shared record engine: gauges x chronological years, no detrending.
rb_matrix <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)) || any(dim(x) == 0L)) {
    stop("Input must be a non-empty numeric matrix or data frame.")
  }
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop("All input columns must be numeric.")
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || any(is.infinite(x))) {
    stop("Input must contain numeric, finite observations or NA.")
  }
  x
}

rb_year <- function(year, n, regular = FALSE) {
  if (!is.numeric(year) || length(year) != n || any(!is.finite(year)) ||
      any(diff(year) <= 0)) {
    stop("Length of year must match columns; years must be finite and strictly increasing.")
  }
  if (regular && length(year) > 2L &&
      !isTRUE(all.equal(diff(year), rep(diff(year)[1], length(year) - 1L)))) {
    stop("SSA requires regularly spaced years; supply missing years explicitly.")
  }
  invisible(year)
}

rb_prepare <- function(list_df, year) {
  if (!is.list(list_df) || is.data.frame(list_df) || !length(list_df)) {
    stop("list_df must be a non-empty list of matrices or data frames.")
  }
  out <- lapply(list_df, rb_matrix)
  if (!all(vapply(out, function(x) identical(dim(x), dim(out[[1]])), logical(1)))) {
    stop("Elements in list_df have different shapes.")
  }
  rb_year(year, ncol(out[[1]]))
  out
}

rb_positive_integer <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(label, " must be a positive integer.")
  }
}

# O(gauges * years). Index input directly instead of copying shuffled matrices.
# Compare current values before updating running maxima and valid counts.
rb_scan <- function(x, perm = seq_len(ncol(x)), detailed = FALSE) {
  p <- nrow(x)
  n <- ncol(x)
  past_max <- rep(-Inf, p)
  past_n <- integer(p)
  observed <- expected <- numeric(n)
  eligible <- integer(n)
  if (detailed) Robs <- R1n <- matrix(NA_real_, p, n, dimnames = dimnames(x))
  for (j in seq_len(n)) {
    now <- x[, perm[j]]
    valid <- !is.na(now)
    use <- valid & past_n > 0L
    obs <- as.numeric(now[use] > past_max[use])
    exp <- 1 / (past_n[use] + 1)
    observed[j] <- sum(obs)
    expected[j] <- sum(exp)
    eligible[j] <- sum(use)
    if (detailed) {
      Robs[use, j] <- obs
      R1n[use, j] <- exp
    }
    past_max[valid] <- pmax(past_max[valid], now[valid])
    past_n[valid] <- past_n[valid] + 1L
  }
  if (detailed) return(list(Robs = Robs, R1n = R1n))
  list(Robs = observed, R1n = expected, n_eligible = eligible)
}

rb_anomaly <- function(observed, expected) {
  out <- rep(NA_real_, length(expected))
  use <- expected > 0
  out[use] <- 100 * (observed[use] - expected[use]) / expected[use]
  out
}

rb_aggregate <- function(matrices, perm = seq_len(ncol(matrices[[1]]))) {
  n <- ncol(matrices[[1]])
  observed <- expected <- numeric(n)
  eligible <- integer(n)
  for (x in matrices) {
    out <- rb_scan(x, perm)
    observed <- observed + out$Robs
    expected <- expected + out$R1n
    eligible <- eligible + out$n_eligible
  }
  data.frame(Robs = observed, R1n = expected,
             Ranom = rb_anomaly(observed, expected), n_eligible = eligible)
}

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
#' Func_RB_matrix(rbind(c(1, 2, 2, 4), c(2, NA, 3, 1)))
Func_RB_matrix <- function(df) rb_scan(rb_matrix(df), detailed = TRUE)

rb_vector <- function(x) {
  if (!is.numeric(x) || !is.null(dim(x)) || any(is.infinite(x))) {
    stop("x must be a numeric vector containing finite observations or NA.")
  }
  if (sum(!is.na(x)) < 2L) stop("At least two valid observations are required.")
  Func_RB_matrix(matrix(x, nrow = 1))
}

#' Identify Record-Breaking Observations
#' @param x Numeric vector with at least two valid observations, in time order.
#' @return Numeric vector: 1 for strict records, 0 for non-records, NA for
#'   missing observations and the first valid observation.
#' @export
Func_Ranom <- function(x) as.numeric(rb_vector(x)$Robs)

#' Expected Record Probabilities
#' @inheritParams Func_Ranom
#' @return Numeric vector of 1/n probabilities; missing observations and the
#'   first valid observation return NA. Missing values do not advance n.
#' @export
Func_Ranom_1N <- function(x) as.numeric(rb_vector(x)$R1n)

#' Aggregate Record Counts Within a Region
#' @inheritParams Func_RB_matrix
#' @param year Numeric, finite, strictly increasing years, one per column.
#' @param plot Logical; print the anomaly plot when TRUE.
#' @return Data frame with year, R1n, Robs, Ranom (percent), and n_eligible.
#'   Ineligible years have zero counts and NA anomaly. Calculations retain
#'   full precision. An eligible observation is valid and has prior history.
#' @export
#' @examples
#' Func_RB_region(rbind(1:4, 4:1), 2001:2004)
Func_RB_region <- function(df, year, plot = FALSE) {
  matrices <- rb_prepare(list(df), year)
  out <- rb_aggregate(matrices)
  out <- data.frame(year = year, out[c("R1n", "Robs", "Ranom", "n_eligible")])
  if (plot) plot_Ranom(out$Ranom, year)
  out
}

#' Aggregate Counts Across Aligned Datasets
#' @param list_df Non-empty list of numeric gauge-by-year matrices or data
#'   frames of identical dimensions. Elements represent aligned months or
#'   comparable series within ONE regional window, not independent regions.
#' @inheritParams Func_RB_region
#' @return List of per-dataset regional results and out_Ranom, containing
#'   time, Robs, R1n, Ranom (percent), and n_eligible. Counts are pooled before
#'   normalization; n_eligible counts gauge-series contributions.
#' @export
Func_RB_region_list <- function(list_df, year, plot = FALSE) {
  matrices <- rb_prepare(list_df, year)
  out <- lapply(matrices, function(x) Func_RB_region(x, year))
  if ("out_Ranom" %in% names(out)) stop("out_Ranom is a reserved list name.")
  out$out_Ranom <- data.frame(time = year, rb_aggregate(matrices))
  if (plot) plot_Ranom(out$out_Ranom$Ranom, year)
  out
}

# Historical internal entry points also use the shared engine.
Func_prepare_matrix_list <- function(list_df) {
  rb_prepare(list_df, seq_len(ncol(list_df[[1]])))
}
Func_RB_matrix_sums <- function(x, perm = NULL) {
  if (is.null(perm)) perm <- seq_len(ncol(x))
  rb_scan(x, perm)
}
Func_RB_region_list_Ranom_only <- function(list_mat, perm = NULL) {
  if (is.null(perm)) perm <- seq_len(ncol(list_mat[[1]]))
  rb_aggregate(list_mat, perm)
}

rb_statistic <- function(x, year, statistic, L, groups, span) {
  switch(statistic,
    raw = x,
    ssa = ssa_trend(x, L = L, groups = groups),
    loess = calculate_loess_trend(data.frame(time = year, value = x), span)
  )
}

rb_quantiles <- function(simulations) {
  q <- vapply(seq_len(ncol(simulations)), function(j) {
    x <- simulations[, j]
    if (all(is.na(x))) return(rep(NA_real_, 4))
    stats::quantile(x, c(0.05, 0.95, 0.025, 0.975), na.rm = TRUE, names = FALSE)
  }, numeric(4))
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
#' @inheritParams Func_RB_region_list
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
  statistic <- match.arg(statistic)
  matrices <- rb_prepare(list_df, year)
  rb_positive_integer(runs, "runs")
  rb_year(year, ncol(matrices[[1]]), regular = statistic == "ssa")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      stop("seed must be a nonnegative integer or NULL.")
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  observed <- rb_aggregate(matrices)
  observed_statistic <- rb_statistic(observed$Ranom, year, statistic, L, groups, span)
  n <- length(year)
  simulated <- matrix(NA_real_, runs, n)
  for (r in seq_len(runs)) {
    perm <- sample.int(n)
    values <- rb_aggregate(matrices, perm)$Ranom
    simulated[r, ] <- tryCatch(
      rb_statistic(values, year, statistic, L, groups, span),
      error = function(e) stop("Simulation ", r, ": ", conditionMessage(e), call. = FALSE)
    )
  }
  list(
    data = data.frame(time = year, observed, statistic_value = observed_statistic,
                      rb_quantiles(simulated)),
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

rb_legacy_ci <- function(fit, loess = FALSE) {
  x <- fit$data
  trend <- if (loess) calculate_loess_trend(data.frame(time = x$time, value = x$Ranom))
    else x$statistic_value
  data.frame(time = x$time, value = x$Ranom, loess_trend = trend,
    x[c("ci_lower_90", "ci_upper_90", "ci_lower_95", "ci_upper_95")])
}

#' Legacy Raw-Anomaly Bands for Aligned Datasets
#' @inheritParams rb_anomaly_ci
#' @return Data frame with time, value, loess_trend, and 90/95 percent bands.
#'   Bands apply to raw anomalies, NOT the accompanying observed LOESS trend.
#' @export
Func_RB_region_list_CI <- function(list_df, year, seed = 1, runs = 1000) {
  rb_legacy_ci(rb_anomaly_ci(list_df, year, "raw", runs, seed), loess = TRUE)
}

#' Legacy Raw-Anomaly Bands for One Dataset
#' @inheritParams Func_RB_matrix
#' @inheritParams rb_anomaly_ci
#' @return Data frame as in [Func_RB_region_list_CI()]. Uses the caller's RNG.
#' @export
RB_Ranom_CI <- function(df, year, runs = 1000) {
  Func_RB_region_list_CI(list(df), year, seed = NULL, runs = runs)
}

#' Legacy SSA-Trend Null Bands
#' @inheritParams rb_anomaly_ci
#' @return Data frame with time, value, loess_trend, and 90/95 percent bands.
#'   For compatibility loess_trend contains SSA, and bands apply to that SSA
#'   statistic. Prefer [rb_anomaly_ci()] for explicit target and settings.
#' @export
Func_RB_region_list_CI_fast <- function(list_df, year, seed = 1, runs = 1000,
                                        L = 15, groups = list(trend = 1)) {
  rb_legacy_ci(rb_anomaly_ci(list_df, year, "ssa", runs, seed, L, groups))
}

#' Legacy LOESS-Trend Null Bands
#' @inheritParams Func_RB_matrix
#' @inheritParams rb_anomaly_ci
#' @return Data frame with time, value, loess_trend, ci_lower and ci_upper
#'   (90 percent LOESS-trend null bands). Uses the caller's RNG.
#' @export
LOESS_Trend_CI <- function(df, year, runs = 1000) {
  out <- rb_legacy_ci(rb_anomaly_ci(list(df), year, "loess", runs, seed = NULL))
  names(out)[4:5] <- c("ci_lower", "ci_upper")
  out[1:5]
}

#' Deprecated Calendar-Reindexed Anomaly Bands
#'
#' Retained only for legacy analyses. Reassigns shuffled sequence-position
#' anomalies to the original observation years. This is NOT the sequence-position
#' null reference for the documented method. Use [RB_Ranom_CI()] or
#' [rb_anomaly_ci()] instead; their results intentionally differ.
#' @inheritParams Func_RB_matrix
#' @inheritParams rb_anomaly_ci
#' @return Legacy data frame with time, value, loess_trend and 90/95 percent bands.
#' @export
RB_Ranom_CI_index <- function(df, year, runs = 1000) {
  .Deprecated("rb_anomaly_ci", msg = paste(
    "RB_Ranom_CI_index is deprecated: calendar-reindexed bands are not",
    "sequence-position null bands. Use rb_anomaly_ci(statistic = 'raw')."))
  matrices <- rb_prepare(list(df), year)
  rb_positive_integer(runs, "runs")
  n <- length(year)
  simulated <- matrix(NA_real_, runs, n)
  for (r in seq_len(runs)) {
    perm <- sample.int(n)
    simulated[r, perm] <- rb_aggregate(matrices, perm)$Ranom
  }
  observed <- rb_aggregate(matrices)
  rb_legacy_ci(list(data = data.frame(time = year, observed,
                                      rb_quantiles(simulated))), loess = TRUE)
}
