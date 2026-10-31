# Repository Guidelines

## Project Structure & Module Organization

This is a script-based R toolkit for record-breaking anomalies in hydrologic and climate time series, not an installable R package.

- `R/Func_Anomaly.R`: original anomaly, trend, and confidence-interval functions.
- `R/RB_anomaly.R`: matrix-based record calculations, regional aggregation, and accelerated Monte Carlo routines.
- `R/SSA_trend.R`: singular spectrum analysis helpers.
- `R/Ranom_plot.R`: plotting helpers.
- `data/mega_floods.RData`: example annual maxima and gauge metadata.
- `test_example.R`: exploratory analysis example, not an automated test suite.
- `references/`: supporting papers; `README.md` describes the scientific framework.

## Build, Test, and Development Commands

Run commands from the repository root with R available on PATH:

```powershell
Rscript -e "install.packages(c('tidyverse', 'trend', 'formattable', 'ggplot2', 'Rssa'))"
Rscript -e "invisible(lapply(list.files('R', pattern='[.]R$', full.names=TRUE), parse))"
Rscript -e "setwd('R'); source('RB_anomaly.R')"
```

These install dependencies, check source syntax, and load the matrix implementation respectively. The working-directory change accommodates its relative `source('./SSA_trend.R')` call. No build system or dependency lockfile is configured.

Before running `Rscript test_example.R`, update its hardcoded workspace and input/source paths to match the current layout and review output destinations. The example contains exploratory sections requiring additional objects. `Ranom_plot.R` currently calls `library(ggplot)`; resolve this dependency typo before using it.

## Coding Style & Naming Conventions

Use UTF-8 and two-space indentation, matching the RStudio project settings. Preserve existing public function names such as `Func_RB_matrix`; use descriptive snake_case names for new helpers and variables. Document matrix orientation, missing-value handling, units, and return fields. No formatter or linter is configured.

## Testing Guidelines

No automated testing framework or coverage threshold exists. Validate changed calculations with small deterministic matrices covering increasing records, ties, missing values, and empty observation years. Preserve gauge-by-year orientation in computational inputs and year alignment. Use fixed seeds and small simulation counts for smoke checks; report simulation settings for substantive validation. Avoid sourcing overlapping implementations together because function names can collide.

## Commit & Pull Request Guidelines

History uses short, topic-focused subjects such as `add plot functions` and `SSA trend estimation`; no strict convention is evident. Keep commits focused. PRs should describe changed behavior, affected functions, validation commands, and scientific assumptions. Include plots when visualization changes. Preserve unrelated local edits and avoid committing generated outputs or credentials.
