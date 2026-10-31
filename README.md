# Record-Breaking (RB) Anomaly

An R toolkit to compute and analyze record-breaking (RB) anomalies in climate time series, with a focus on hydrologic extremes.

## Local R package

The installable package is named **RBanomaly** (version 0.1.0). It is intended
for local use; no CRAN publication is needed. The repository directory can
remain named `RB_anomaly`.

Install the runtime dependencies once from an R session:

```r
install.packages(c("ggplot2", "trend", "formattable", "Rssa"))
```

From a terminal in this repository, install the package:

```powershell
R.exe CMD INSTALL .
```

Use `R CMD INSTALL .` on Linux/macOS. The selected library must be writable;
use `--library=path/to/library` to choose one explicitly. Alternatively, in R:

```r
install.packages("remotes") # once, if needed
remotes::install_local(".", dependencies = NA, upgrade = "never")
```

Start an analysis in a fresh R session:

```r
library(RBanomaly)
x <- rbind(c(1, 2, 2, 4), c(2, NA, 3, 1)) # stations x years
result <- Func_RB_region(x, year = 2001:2004)
plot_Ranom(result$Ranom, result$year)
?Func_RB_region
help(package = "RBanomaly")
```

The first valid observation is excluded from record comparisons; ties are not
new records. Missing values do not increase the valid observation count.
Always supply years in the same order as matrix columns.

Load the existing example data through the package:

```r
data("mega_floods", package = "RBanomaly")
year <- df_AMS_All$year
x <- t(as.matrix(df_AMS_All[, 2:4])) # transpose gauge columns
result <- Func_RB_region(x, year)
```

## Unified record and simulation workflow (31 October 2026)

All record calculations now live in `R/RB_anomaly.R`. Vector, matrix,
regional, and Monte Carlo functions share one running-maximum/count engine.
The simulation path aggregates directly without allocating gauge-level output
matrices or copying entire shuffled datasets.

Use the explicit simulation interface for new analyses:

```r
library(RBanomaly)
set.seed(42)
x <- matrix(rexp(4 * 40), nrow = 4) # four gauges x forty years
fit <- rb_anomaly_ci(
  list(x), year = 1981:2020,
  statistic = "ssa", runs = 1000, seed = 1,
  L = 15, groups = list(trend = 1)
)
fit$data       # counts, eligibility, anomaly, selected statistic, null bands
fit$settings   # statistic, seed, simulation and smoothing settings
plot_Ranom_CI(fit)
```

`statistic = "ssa"` (default), `"raw"`, or `"loess"` selects the quantity
compared with its null distribution. Apply SSA/LOESS to **each simulation
before calculating quantiles**. `data$statistic_value` and the interval bounds
refer to the same quantity; `data$Ranom` always contains the raw anomaly.
The 90% and 95% bands are pointwise null-reference ranges, not confidence
intervals around the observed trend and not simultaneous significance bands.
For a raw analysis, no smoother or minimum smoothing length is required.

Each simulation permutes year columns synchronously across all gauges and
all aligned list elements. A list represents months or comparable series in
**one regional window**, not independently windowed geographic regions.
Missing entries move with their columns; valid history and expectations are
recomputed for every permutation. Temporal dependence is destroyed. This
release does not select regional windows, test iid, or perform block shuffling.
These assumptions must be assessed for a discharge application.

The paper uses 10,000 simulations; the default remains 1,000 for compatibility.
Set `runs` explicitly for substantive work and assess Monte Carlo precision.
SSA's `L = 15` follows the paper's window, but `groups = list(trend = 1)` is a
package choice, not a verified author component-selection rule. SSA requires
regularly spaced years and sufficient usable length. Leading/trailing missing
values remain missing; interior gaps are rejected, including gaps arising in a
simulation (the error identifies its run). No interpolation is performed.
LOESS uses `span = 0.8` by default, requires a usable smoothing neighborhood,
and preserves missing positions. No observed discharge series is detrended
before counting its records.

### Calculation and compatibility changes

- Strict records, exclusion of the first valid observation, and the valid-count
  `1/n` expectation are unchanged. Pool counts before normalizing; the first
  observation is excluded from both numerator and denominator.
- All calculation paths now retain full precision. Earlier regional functions
  rounded anomalies to one decimal. Round explicitly only for display.
- Regional outputs add `n_eligible`: valid gauge-series observations with prior
  history. Years without eligible observations have `Robs = R1n = 0` and
  `Ranom = NA`; the fast path previously returned missing counts.
- Numeric inputs are validated without coercing text or factors into numbers.
  Infinite observations, empty inputs, incompatible shapes, and invalid years
  or simulation settings produce errors. Matrix dimnames are retained.
- Existing public function names and positional arguments remain available.
  Vector functions still reject fewer than two valid observations; matrices
  accommodate such gauges as ineligible contributors.
- Seeded simulations are reproducible and restore the caller's random-number
  state. This changes the previous seed side effect. `seed = NULL` consumes
  the caller's RNG stream; wrappers without a seed argument retain that behavior.
- The exploratory example loads `RBanomaly` rather than sourcing removed files.
  Install/reload the package after updating; direct sourcing of
  `Func_Anomaly.R` is no longer supported.

| Existing interface | Preserved statistical target and output |
| --- | --- |
| `Func_Ranom()`, `Func_Ranom_1N()` | Vector record indicators / expected probabilities |
| `Func_RB_matrix()` | Gauge-level record indicators and expectations |
| `Func_RB_region()`, `Func_RB_region_list()` | Pooled counts and full-precision anomalies; adds eligibility |
| `RB_Ranom_CI()`, `Func_RB_region_list_CI()` | Raw-anomaly bands; accompanying `loess_trend` is descriptive, not the band target |
| `Func_RB_region_list_CI_fast()` | SSA-trend bands; legacy `loess_trend` field still contains SSA |
| `LOESS_Trend_CI()` | LOESS-trend bands, with legacy 90% `ci_lower` / `ci_upper` fields |
| `RB_Ranom_CI_index()` | Deprecated calendar-reindexed calculation; retained with a warning, excluded from the recommended workflow |

`plot_Ranom_CI()` accepts the new result, its `data` component, and legacy
outputs (including 90%-only LOESS results). For new results it draws the
selected statistic against matching bands. With legacy raw-band results,
do not interpret their separate observed LOESS line against the raw bands.

### Source organization and validation

- `R/RB_anomaly.R`: shared record engine, aggregation, simulation, and compatibility wrappers.
- `R/SSA_trend.R`: generic SSA reconstruction helpers.
- `R/trend_utils.R`: Sen slope, linear slope, Mann-Kendall, and LOESS helpers.
- `R/Ranom_plot.R`: plots and table formatters.
- `R/package.R`, `R/data.R`: package imports and example-data documentation.

`R/Func_Anomaly.R` has been removed after relocating its functions. The public
exports remain available, with `rb_anomaly_ci()` added. Tests cover exact
record probabilities by permutation enumeration, ties and missing histories,
empty eligible years, count pooling, full/aggregate equivalence, manual
raw/SSA/LOESS null calculations, seeded reproducibility, compatibility, and
plot construction. See the development commands below to rerun validation.

The [methodology guide](references/record_breaking_anomaly_methodology.md)
distinguishes the paper's methods from package choices and discharge-specific
assumptions. Consolidation does not establish iid behavior in a new dataset or
constitute complete reproduction of the paper's multi-region analysis.

## Package development

- `DESCRIPTION` declares metadata and dependencies; `NAMESPACE` declares exports.
- `R/` contains functions and roxygen2 documentation; `man/` contains generated help.
- `data/` contains the existing megaflood example; `tests/testthat/` contains regression tests.
- `references/` and the exploratory `test_example.R` remain repository resources
  and are excluded from the built package. The exploratory script still uses
  machine-specific paths; use the examples above for the package workflow.

From the repository root:

```r
install.packages(c("roxygen2", "testthat")) # development only
roxygen2::roxygenise(".")                  # regenerate man/ and NAMESPACE
testthat::test_local(".")                 # run regression tests
```

To build and check in a terminal (use `R` instead of `R.exe` on Linux/macOS):

```powershell
R.exe CMD build --resave-data .
R.exe CMD check --no-manual RBanomaly_0.1.0.tar.gz
```

`--resave-data` compresses the data in the built archive without changing the
original dataset. `--no-manual` avoids requiring LaTeX; help pages and their
examples are still checked. Edit roxygen comments in `R/`, then regenerate documentation. Use
the installed namespace instead of sourcing individual files. No open-source
redistribution license was selected; see `LICENSE` for the local-use notice.

This repository implements methods for detecting changes in the frequency of record-breaking events relative to a stationary climate. 


# Overview

Extreme climate events, especially record-breaking precipitation and temperature, have increased in frequency under global warming. Instead of relying on parametric extreme value theory, the record-breaking framework offers a robust, distribution-free approach.

A value is considered record-breaking if it exceeds all previous values in a time series. The key metric implemented here is the Record-Breaking Anomaly (RB anomaly):

$$
RB_{\text{anom}} = \frac{
    R_{\text{obs}} - R_{\text{exp}}
}{R_{\text{exp}} } \times 100 \%
$$

where:

- $R_{\text{obs}}$: observed number of record-breaking events
- For one complete series, $R_{\text{exp}} = \sum_{n=2}^{N}{\frac{1}{n}}$:
  expected count under the continuous iid reference, excluding the first
  observation. For a regional year or reporting period, sum the eligible
  observations' individual $1/n$ probabilities over exactly the same support
  as the observed count; $n$ is each series' current valid-observation count.

This method allows direct comparison between observations and a null hypothesis of no climate trend.

# Features

- Compute record-breaking events for time series (e.g., precipitation, temperature)
- Estimate RB anomaly relative to iid expectation
- Perform Monte Carlo shuffling to generate confidence intervals
- Spatial aggregation across regions


# Reference
1. Lehmann, J., Coumou, D. and Frieler, K.  2015.  Increased record-breaking precipitation events under global warming. Climatic Change, 132(4), 501-515. doi: 10.1007/s10584-015-1434-y.
2. Robinson, A., Lehmann, J., Barriopedro, D., Rahmstorf, S. and Coumou, D.  2021.  Increasing heat and rainfall extremes now far outside the historical climate. npj Climate and Atmospheric Science, 4(1), 45. doi: 10.1038/s41612-021-00202-w.

# Contact

[Xiaoxiang Guan](guan@gfz.de)

[GFZ Helmholtz Centre for Geosciences](https://www.gfz.de)

