#' Record-Breaking Anomaly Analysis
#'
#' Functions for record indicators, regional anomalies, temporal-shuffling
#' confidence intervals, and trend visualization. Computational matrices have
#' stations in rows and ordered observation years in columns.
#'
#' @importFrom stats coefficients lm loess predict quantile
#' @importFrom trend sens.slope mk.test
#' @importFrom formattable formatter style
#' @importFrom utils globalVariables
#' @importFrom ggplot2 ggplot aes geom_col geom_hline geom_ribbon geom_line labs theme_minimal
#' @keywords internal
"_PACKAGE"

# Names evaluated within ggplot2 aesthetics.
utils::globalVariables(c("time", "value", "year", "Ranom", "loess_trend",
                        "ci_lower_90", "ci_upper_90", "ci_lower_95", "ci_upper_95"))
