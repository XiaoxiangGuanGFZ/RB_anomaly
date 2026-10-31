#' Estimate Sen Slope
#' 
#' Use trend::sens.slope() with a 95 percent confidence level.
#' 
#' @param x Numeric time series.
#' @return Numeric vector containing the estimated slope and p-value.
#' @export
senslope <- function(x) 
{
  # ---- description
  # a function to test the trend with Sen's slope method
  # ---- input:
  # x: a multiple value vector, a time series
  # ---- output
  # slope and the p.value from test
  out = sens.slope(x, conf.level = 0.95)
  
  return(
    c(out$estimates, out$p.value)
  )
}

#' Estimate a Linear Trend
#' 
#' Fit a linear regression against the observation index.
#' 
#' @param x Numeric time series.
#' @return Named numeric scalar containing the regression slope.
#' @export
lmslope <- function(x) {
  # ---- description:
  # linear regression method-derived trend
  # 
  # ---- input:
  # x: a multiple value vector, a time series
  # ---- output
  # slope
  # 
  t = seq(1, length(x), 1)
  df = data.frame(
    t, x
  )
  reg = lm(x~t, data = df)
  slopecoeff = coefficients(reg)
  return(slopecoeff[2])
}

#' Test for a Monotonic Trend
#' 
#' Perform the Mann-Kendall test using trend::mk.test().
#' 
#' @param x Numeric time series.
#' @return Numeric vector containing the test statistic and p-value.
#' @export
mk_test <- function(x)
{
  # ---- description:
  # trend detection with Mann-Kendall rank method
  # 
  # ---- input:
  # x: a multiple value vector, a time series
  # ---- output
  # MK statistic and p.value
  # 
  mk_out = mk.test(x)
  mk_z = mk_out$statistic[1]
  return(
    c(mk_z, mk_out$p.value)
  )
}

#' Smooth an Anomaly Series with LOESS
#' 
#' Fit value against time and predict at the original times.
#' 
#' @param data Data frame containing numeric time and value columns.
#' @param span LOESS smoothing span.
#' @return Numeric vector of fitted values at the input times.
#' @export
calculate_loess_trend <- function(data, span = 0.8) {
  if (!is.data.frame(data) || !all(c("time", "value") %in% names(data)) ||
      !is.numeric(data$time) || !is.numeric(data$value) ||
      any(!is.finite(data$time)) || any(is.infinite(data$value))) {
    stop("data must contain finite numeric time and numeric value (NA allowed).")
  }
  if (!is.numeric(span) || length(span) != 1L || !is.finite(span) || span <= 0) {
    stop("span must be a positive finite number.")
  }
  out <- rep(NA_real_, nrow(data))
  valid <- which(!is.na(data$value))
  if (!length(valid)) return(out)
  if (length(unique(data$time[valid])) < 4L || floor(length(valid) * span) < 4L) {
    stop("LOESS requires at least four distinct valid times in its smoothing neighborhood.")
  }
  loess_model <- stats::loess(value ~ time, data = data, span = span)
  out[valid] <- stats::predict(loess_model, data.frame(time = data$time[valid]))
  if (any(!is.finite(out[valid]))) stop("LOESS produced non-finite fitted values.")
  out
}

