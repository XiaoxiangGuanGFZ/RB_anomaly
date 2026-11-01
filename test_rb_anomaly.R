# ============================================================================
# RBanomaly walkthrough: annual maximum discharge -> records -> null bands
# ============================================================================
# This is a runnable teaching example, not a testthat test or package source.
# Work through the numbered sections in RStudio, or run from the repository:
#   R CMD INSTALL .                      # install the current checkout first
#   Rscript test_rb_anomaly.R
# In an R session after installation:
#   source("test_rb_anomaly.R", echo = FALSE)
# Restart R/reload the package after reinstalling to avoid an older namespace.
# Runtime dependencies and installation instructions are in README.md.
# No packages are installed automatically and no files are written by default.
# Read references/record_breaking_anomaly_methodology.md alongside this script.
# Main examples use the unified rb_* API. Older public names remain aliases;
# section 9 alone demonstrates the older CI interfaces for existing analyses.

# ---- 1. Choose demonstration settings ---------------------------------------
# These are illustrative choices, not thresholds prescribed by the paper.
# All settings that affect selection, simulation, and saving are collected here.
config <- list(
  start_year = 1961,
  end_year = 2010,
  region_code = 1,             # df_gauge$regions1 code; no inferred region name
  max_area_km2 = 1000,
  dams = "no",                # use the dataset's recorded category
  max_gauges = 20L,            # deterministic subset for a quick walkthrough
  runs = 100L,                 # demonstration only, not publication precision
  seed = 42L,
  ssa_window = 15L,
  ssa_groups = list(trend = 1), # package choice; not a verified author recipe
  loess_span = 0.8,
  run_legacy_examples = TRUE,
  legacy_runs = 20L,           # small runs only to illustrate old interfaces
  show_plots = interactive(),  # Rscript will not silently create Rplots.pdf
  save_outputs = FALSE,
  output_dir = file.path("outputs", "rb_walkthrough"),
  overwrite = FALSE
)

if (!requireNamespace("RBanomaly", quietly = TRUE)) {
  stop("Install this checkout first with R CMD INSTALL . (see README.md).")
}
library(RBanomaly)
if (!all(c("rb_anomaly_ci", "rb_record_matrix", "rb_plot_anomaly_ci") %in%
         getNamespaceExports("RBanomaly"))) {
  stop("Installed RBanomaly is outdated. Reinstall this checkout and restart R.")
}
if (config$start_year >= config$end_year || config$max_gauges < 1L ||
    config$max_gauges != floor(config$max_gauges)) {
  stop("Choose an increasing year window and a positive integer max_gauges.")
}
message("[1/10] Package: ", packageVersion("RBanomaly"),
        " from ", find.package("RBanomaly"))
print(config)

# ---- 2. Load data and align gauge metadata ----------------------------------
# Use a separate environment so loading example objects does not overwrite
# similarly named objects in your interactive workspace.
example_data <- new.env()
data("mega_floods", package = "RBanomaly", envir = example_data)
ams <- example_data$df_AMS_All
gauges <- example_data$df_gauge

# ams: rows = years, first column = year, remaining columns = gauge series.
# gauges: id, lat, lon, area, regions1, regions2, dams, plus original index X.
# Values are annual maximum SPECIFIC discharge, documented as m3/s/km2.
# The RB anomaly is dimensionless (%): it measures record frequency, not the
# size of a flood or the amount by which the previous record was exceeded.
message("[2/10] Bundled data: ", nrow(ams), " years, ", ncol(ams) - 1L,
        " gauges; calendar ", min(ams$year), "--", max(ams$year))
print(head(gauges))

# Explicit identity mapping: metadata id 123 corresponds to AMS column s123.
# Never assume that metadata row number is the gauge's input matrix row.
gauges$series_name <- paste0("s", as.character(gauges$id))
stopifnot(!anyDuplicated(gauges$series_name), !anyDuplicated(names(ams)),
          !anyDuplicated(ams$year), all(diff(ams$year) > 0),
          all(gauges$series_name %in% names(ams)))

# ---- 3. Define the sample and convert to gauges x years ---------------------
year_rows <- which(ams$year >= config$start_year & ams$year <= config$end_year)
year <- ams$year[year_rows]
if (!identical(as.numeric(year), as.numeric(seq(config$start_year, config$end_year)))) {
  stop("Requested window must be present as consecutive annual observations.")
}
candidate_rows <- which(gauges$regions1 == config$region_code &
                          gauges$area <= config$max_area_km2 &
                          gauges$dams == config$dams)
candidate_metadata <- gauges[candidate_rows, , drop = FALSE]
candidate_values <- ams[year_rows, candidate_metadata$series_name, drop = FALSE]

# For this teaching example retain only gauges with finite values in EVERY
# selected year. This makes SSA valid for every permutation: apart from the
# excluded first year, the aggregated anomaly has no missing interior years.
# This is a sample-selection choice, NOT a package requirement for raw counts.
# Complete-case selection can bias coverage; assess it in a real study.
complete <- vapply(candidate_values, function(z) all(is.finite(z)), logical(1))
selected_metadata <- candidate_metadata[complete, , drop = FALSE]
if (!nrow(selected_metadata)) {
  stop("No complete gauges for these settings. Review the region/year window.")
}
# Select deterministically by gauge identifier; do not rank by anomaly magnitude.
selected_metadata <- selected_metadata[order(selected_metadata$id), , drop = FALSE]
selected_metadata <- head(selected_metadata, config$max_gauges)
x <- t(as.matrix(ams[year_rows, selected_metadata$series_name, drop = FALSE]))
colnames(x) <- as.character(year)
stopifnot(identical(rownames(x), selected_metadata$series_name),
          ncol(x) == length(year), all(is.finite(x)))
if (config$ssa_window >= length(year) - 1L) {
  stop("SSA window must be smaller than the usable anomaly length (years minus one).")
}
message("[3/10] Candidates: ", nrow(candidate_metadata),
        "; complete: ", sum(complete), "; selected: ", nrow(x),
        "; years: ", min(year), "--", max(year))
print(selected_metadata[, c("id", "series_name", "area", "regions1", "dams")])

# IMPORTANT: slicing BEFORE record counting makes start_year the beginning of
# record history. These are records relative to that window, not necessarily
# records relative to all observations since 1805. To use longer history,
# calculate on that longer matrix first and subset the results for reporting.
# This small sample is not a claim about all gauges in the metadata region.

# ---- 4. Inspect records at one gauge and across the matrix ------------------
one_gauge <- x[1, ]
record_indicator <- rb_record_indicator(one_gauge)
record_probability <- rb_record_probability(one_gauge)
gauge_example <- data.frame(
  year = year,
  specific_discharge = as.numeric(one_gauge),
  Robs = record_indicator,
  R1n = record_probability
)
message("[4/10] First gauge: ", rownames(x)[1], "; first ten comparisons")
print(head(gauge_example, 10))

# Robs = 1: strict new record; 0: valid observation that did not set a record.
# R1n = 1/n: iid expected probability at valid observation number n, not a
# return-period estimate. First valid observations have NA in both columns.
# For this complete series: NA, 1/2, 1/3, 1/4, ... . Ties are not new records.
records <- rb_record_matrix(x)
stopifnot(isTRUE(all.equal(as.numeric(records$Robs[1, ]), record_indicator)),
          isTRUE(all.equal(as.numeric(records$R1n[1, ]), record_probability)))

# Missing-data demonstration, derived from this same observed gauge:
# remove its third value to see that missing years do not advance n. Removing
# an observation can also change later records if it was the previous maximum.
with_gap <- one_gauge
with_gap[3] <- NA_real_
missing_example <- data.frame(
  year = head(year, 6),
  discharge = head(with_gap, 6),
  Robs = head(rb_record_indicator(with_gap), 6),
  R1n = head(rb_record_probability(with_gap), 6)
)
print(missing_example)
# Expected probabilities now begin NA, 1/2, NA, 1/3, 1/4, 1/5.
# We do not use this artificially modified series for the main analysis.

# ---- 5. Pool counts and interpret the regional anomaly ----------------------
regional <- rb_anomaly_region(x, year, plot = FALSE)
message("[5/10] Regional counts and anomaly")
print(head(regional, 10))
# R1n: sum of expected probabilities across eligible gauges in this year.
# Robs: observed count of gauges setting records in this year.
# n_eligible: gauges with a current value AND at least one earlier valid value.
# Ranom = 100 * (Robs - R1n) / R1n, retaining full numerical precision.
# +50% means 1.5 times the expected record count; -50% means half as many.
# -100% means zero observed records despite a positive expectation.
# With no eligible gauges, counts are zero and Ranom is NA (not 0%).
stopifnot(isTRUE(all.equal(regional$Robs, unname(colSums(records$Robs, na.rm = TRUE)))),
          isTRUE(all.equal(regional$R1n, unname(colSums(records$R1n, na.rm = TRUE)))))

# A list interface supports aligned series (e.g. separate calendar months).
# This dataset has annual maxima only, so use list(x), not invented months.
# A list of independently windowed geographical regions is NOT supported by
# the single shared-permutation workflow. Do not duplicate x to inflate counts.
list_result <- rb_anomaly_region_list(list(annual_maxima = x), year)
stopifnot(isTRUE(all.equal(list_result$out_Ranom$Ranom, regional$Ranom)))
print(head(list_result$out_Ranom, 4))
# $annual_maxima holds its individual result; $out_Ranom pools counts first.

# A pooled period anomaly differs from an unweighted mean of annual anomalies.
period_observed <- sum(regional$Robs)
period_expected <- sum(regional$R1n)
period_summary <- data.frame(
  start_year = min(year), end_year = max(year),
  Robs = period_observed, R1n = period_expected,
  Ranom = 100 * (period_observed - period_expected) / period_expected
)
print(period_summary)
# This summary is descriptive. The pointwise bands below do not test this
# period statistic; that would require computing it for each shuffled run.

# ---- 6. Simulate RAW-anomaly null-reference bands ---------------------------
message("[6/10] Raw null bands: ", config$runs, " runs; seed ", config$seed)
raw_fit <- rb_anomaly_ci(list(x), year, statistic = "raw",
                         runs = config$runs, seed = config$seed)
print(head(raw_fit$data, 6))
print(raw_fit$settings)
# $data includes observed counts, raw Ranom, statistic_value, and 90/95% bounds.
# In raw mode, statistic_value equals Ranom. Each run permutes whole year
# columns together across all gauges, preserving same-year spatial association.
# Record histories and expectations are recomputed in the new sequence order.
# This is an iid/exchangeable-year reference, not proof that discharge is iid.

# ---- 7. Match SSA/LOESS trends to their OWN simulated bands ------------------
message("[7/10] SSA and LOESS null bands: ", config$runs, " runs each")
ssa_fit <- rb_anomaly_ci(
  list(x), year, statistic = "ssa", runs = config$runs, seed = config$seed,
  L = config$ssa_window, groups = config$ssa_groups
)
loess_fit <- rb_anomaly_ci(
  list(x), year, statistic = "loess", runs = config$runs, seed = config$seed,
  span = config$loess_span
)
print(head(ssa_fit$data[, c("time", "Ranom", "statistic_value",
                           "ci_lower_95", "ci_upper_95")], 6))
# In SSA mode, Ranom is still raw, but statistic_value is the observed SSA trend.
# Bands come from SSA applied to EACH shuffled anomaly before taking quantiles.
# Compare the SSA line with SSA bands; do not judge raw bars against those bands.
# The same rule applies to LOESS. These are distinct statistical targets.
# L=15 and the selected eigentriples must be reported and checked for sensitivity.
# SSA removes/restores missing edges, but rejects interior gaps rather than
# inventing values. Both observations and simulations must satisfy this rule.

# Direct smoothers reproduce the observed statistics; they do not alone provide
# null bands. rb_ssa_trend_matrix() applies SSA to each row (series x time).
observed_ssa <- rb_ssa_trend(regional$Ranom, L = config$ssa_window,
                           groups = config$ssa_groups)
observed_ssa_matrix <- rb_ssa_trend_matrix(
  rbind(regional_anomaly = regional$Ranom),
  L = config$ssa_window, groups = config$ssa_groups
)
observed_loess <- rb_loess_trend(
  data.frame(time = year, value = regional$Ranom), span = config$loess_span
)
stopifnot(isTRUE(all.equal(observed_ssa, ssa_fit$data$statistic_value)),
          isTRUE(all.equal(as.numeric(observed_ssa_matrix[1, ]), observed_ssa)),
          isTRUE(all.equal(observed_loess, loess_fit$data$statistic_value)))

# ---- 8. Read the comparison and draw plots ----------------------------------
comparison <- ssa_fit$data
comparison$position_95 <- ifelse(
  is.na(comparison$statistic_value) | is.na(comparison$ci_lower_95) |
    is.na(comparison$ci_upper_95),
  "undefined",
  ifelse(comparison$statistic_value > comparison$ci_upper_95, "above",
         ifelse(comparison$statistic_value < comparison$ci_lower_95,
                "below", "within"))
)
message("[8/10] Observed SSA statistic versus pointwise 95% null bands")
print(table(comparison$position_95))
print(head(comparison[, c("time", "statistic_value", "ci_lower_95",
                          "ci_upper_95", "position_95")], 8))
# "above" / "below": unusually many / few records in the SMOOTHED statistic
# under this specified null. "within" does not prove stationarity or no change.
# These are pointwise bands, not a simultaneous test over all years; neighboring
# smoothed values are dependent. Many year-by-year checks raise multiplicity.
# A departure does not attribute change to climate: regulation, land use,
# measurement changes, and unsuitable dependence assumptions also need review.
# 100 runs makes this an execution demonstration, not a publication conclusion.
# The paper used 10,000; higher runs reduce simulation error, not model error.

plots <- list(
  raw_anomaly = rb_plot_anomaly(regional$Ranom, year, p_print = FALSE),
  raw_bands = rb_plot_anomaly_ci(raw_fit, p_print = FALSE),
  ssa_bands = rb_plot_anomaly_ci(ssa_fit, p_print = FALSE),
  loess_bands = rb_plot_anomaly_ci(loess_fit, p_print = FALSE)
)
if (config$show_plots) for (p in plots) print(p)
# In RStudio you can also run print(plots$ssa_bands) after sourcing this script.
# Internal rb_* helpers implement preparation/scanning/aggregation; applications
# should call the public functions above instead of RBanomaly::: internals.

# ---- 9. Recognize compatibility functions and separate trend diagnostics ----
message("[9/10] Compatibility interfaces and descriptive trend helpers")
if (config$run_legacy_examples) {
  # Two raw-band wrappers: their LOESS column is NOT the interval target.
  legacy_raw_list <- Func_RB_region_list_CI(
    list(x), year, seed = config$seed, runs = config$legacy_runs
  )
  set.seed(config$seed) # wrapper without a seed argument uses caller's RNG
  legacy_raw_single <- RB_Ranom_CI(x, year, runs = config$legacy_runs)
  stopifnot(isTRUE(all.equal(legacy_raw_list, legacy_raw_single)))

  # SSA wrapper: historical field loess_trend actually contains the SSA trend.
  legacy_ssa <- Func_RB_region_list_CI_fast(
    list(x), year, seed = config$seed, runs = config$legacy_runs,
    L = config$ssa_window, groups = config$ssa_groups
  )
  set.seed(config$seed)
  legacy_loess <- LOESS_Trend_CI(x, year, runs = config$legacy_runs)
  print(head(legacy_ssa, 3))
  print(names(legacy_loess)) # ci_lower / ci_upper are 90% LOESS bounds only
  # RB_Ranom_CI_index() is deprecated: it reassigns shuffled-position anomalies
  # to original years. Do not use it for the recommended null comparison.
}

# These generic helpers describe the FIRST gauge's discharge trend, not the
# record anomaly and not the significance of SSA. Here years are consecutive,
# so one observation step is one year and slopes have units (m3/s/km2)/year.
# Removing missing years would break that interpretation for index-based slopes.
linear_slope <- rb_linear_slope(as.numeric(one_gauge))
sen_result <- rb_sen_slope(as.numeric(one_gauge)) # slope and p-value
mk_result <- rb_mann_kendall(as.numeric(one_gauge))  # Mann-Kendall statistic and p-value
print(list(linear_slope = linear_slope, sen_slope_and_p = sen_result,
           mann_kendall_statistic_and_p = mk_result))
# Trend-test p-values need their own dependence assumptions; these calls are
# not iid validation. rb_format_trend()/rb_format_pvalue() only format HTML
# table cells; they do not alter estimates or perform additional testing.

# ---- 10. Optional explicit export and reproducibility record ---------------
# Keep selected IDs and settings with results so the analysis can be reproduced.
walkthrough_results <- list(
  config = config, selected_gauges = selected_metadata, year = year,
  regional = regional, period_summary = period_summary,
  raw = raw_fit, ssa = ssa_fit, loess = loess_fit,
  ssa_comparison = comparison, session_info = sessionInfo()
)
if (config$save_outputs) {
  # Check every destination before writing; overwriting requires an explicit opt-in.
  destinations <- file.path(config$output_dir,
                            c("walkthrough_results.rds", "ssa_comparison.csv"))
  if (!config$overwrite && any(file.exists(destinations))) {
    stop("Output files already exist; choose another output_dir or set overwrite = TRUE.")
  }
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(walkthrough_results, destinations[1])
  write.csv(comparison, destinations[2], row.names = FALSE)
  message("Saved: ", paste(destinations, collapse = ", "))
}
message("[10/10] Walkthrough complete. Results: walkthrough_results; plots: plots.")
message("Files saved: ", config$save_outputs,
        ". Review sample selection and assumptions before scientific interpretation.")

