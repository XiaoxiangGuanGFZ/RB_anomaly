#' Example Annual Maximum Discharges and Gauge Metadata
#'
#' Existing megaflood example data bundled with this repository. The annual
#' maxima table uses a different orientation from computational inputs:
#' transpose its gauge columns before calling the anomaly functions.
#' @name mega_floods
#' @aliases df_AMS_All df_gauge
#' @format Two data frames:
#' \describe{
#'   \item{df_AMS_All}{218 rows and 8024 columns: year followed by 8023 gauge
#'     series of annual maximum specific discharge in cubic metres per second
#'     per square kilometre.}
#'   \item{df_gauge}{8023 rows and 8 columns: X (original index), id (gauge
#'     identifier), lat, lon, area, regions1, regions2, and dams.}
#' }
#' @source Megaflood dataset: \url{https://github.com/tuwhydro/megafloods}.
#'   Metadata and units follow the repository's original example script.
#' @examples
#' data("mega_floods", package = "RBanomaly")
#' year <- df_AMS_All$year
#' x <- t(as.matrix(df_AMS_All[, 2:4]))
#' head(Func_RB_region(x, year))
NULL
