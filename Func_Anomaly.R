# ------------------------ script description -----------------
# Record-Breaking (RB) anomaly: 
# How often extreme records are broken compared to what 
# we would expect in a stable (unchanging) climate.
# 
# A R script used to compute the RB anomaly
# 
# Core references:
# Lehmann, J., Coumou, D. and Frieler, K.  2015.  Increased record-breaking 
#   precipitation events under global warming. Climatic Change, 132(4), 501-515. 
#   doi: 10.1007/s10584-015-1434-y.
# Robinson, A., Lehmann, J., Barriopedro, D., Rahmstorf, S. and Coumou, D.  2021.  
#   Increasing heat and rainfall extremes now far outside the historical climate. 
#   npj Climate and Atmospheric Science, 4(1), 45. doi: 10.1038/s41612-021-00202-w.
# 
# Author: 
# Xiaoxiang Guan (guan@gfz.de)
# Last update:
# 2026-04-21
# 




library(tidyverse)
library(trend)
library(formattable)

# ----------------- Functions: trend estimation ---------------
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

lmslope <- function(x) {
  # ---- description:
  # linear regression method-derived trend
  # 
  # ---- input:
  # x: a multiple value vector, a time series
  # ---- output
  # slope
  # 
  t <- seq(1, length(x), 1)
  df = data.frame(
    t, x
  )
  reg <- lm(x~t, data = df)
  slopecoeff <- coefficients(reg)
  return(slopecoeff[2])
}

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
  mk_out <- mk.test(x)
  mk_z = mk_out$statistic[1]
  return(
    c(mk_z, mk_out$p.value)
  )
}

# ----------- table formatter ------------
# specifier for visualization 
trend_formatter <- formatter(
  "span",
  style = x ~ style(
    color = ifelse(x > 0, "red", "green"), # Green for positive, red for negative
    font.weight = "bold"
  ),
  x ~ sprintf("%+.3f", x) # Format to include + or - sign with one decimal
)

pvalue_formatter <- formatter(
  "span",
  style = x ~ style(
    background.color = ifelse(x < 0.05, "gray", "transparent"), # Gray background for p < 0.05
    color = ifelse(x < 0.05, "white", "black") # White text if background is gray for contrast
  ),
  x ~ sprintf("%.3f", x) # Format p-values to three decimal places
)

# --------------- Functions: record-breaking anomaly ------------

RBF <- function(df, year, plot = FALSE) 
{
  # ---- description:
  # compute the record-breaking anomaly based 
  # on AMS across multiple gauges in a region
  # ---- input:
  # df: a dataframe, with ncol years of observation, and nrow gauges. 
  # year: a vector, indicating the observation year, consistent with columns in `df`
  # ---- output:
  # return a dataframe, 
  # 

  dims = dim(df)
  p = dims[1]  # number of stations
  n = dims[2]  # years
  # make a copy 
  Robs <- matrix(NA, p, n)
  R1n  <- matrix(NA, p, n)
  # define the first year (observation)
  Robs[, 1] = NA
  R1n[, 1] = NA
  
  for (i in 1:p) {
    # iterate each gauge
    for (j in 2:n) {
      # iterate the years from second year
      past = df[i, 1:(j - 1)]
      if (is.na(df[i, j]))
      {
        # this position is NA
        Robs[i, j] = NA 
        R1n[i, j] = NA 
      } else if (sum(!is.na(past)) == 0) {
        # this position is not NA, and
        # the records before this position are all NA
        Robs[i, j] = NA 
        R1n[i, j] = NA
      } else {
        # this position is not NA, and 
        # the records before this position are not all NA,
        # sum(!is.na(past)) > 0
        Robs[i, j] = ifelse(df[i, j] > max(past, na.rm = TRUE), 1, 0)
        R1n[i, j] = 1 / (sum(!is.na(past)) + 1)
      }
    }
  }

  # ensure consistent support
  mask = !is.na(Robs)
  R1n[!mask] = NA
  
  out = data.frame(
    year = year,
    R1n  = colSums(R1n, na.rm = TRUE),
    Robs = colSums(Robs, na.rm = TRUE)
  )

  # function output:
  # a data.frame with four columns: year, R1n, Robs, Ranom
  out$Ranom = ifelse(out$R1n > 0,
                     round((out$Robs - out$R1n)/out$R1n * 100, 1),
                     NA)
  
  # visualize the RB anomaly
  if (plot)
  {
    ggplot(out) + geom_bar(aes( x = year, y = Ranom), stat = "identity")
  }

  return(out)
}


RB_Ranom_CI <- function(df, year, runs = 1000){
  # ---- description:
  # Compute confidence intervals (CI) for the record-breaking anomaly (Ranom)
  # using a Monte Carlo null model based on temporal shuffling (iid assumption).
  #
  # The method follows Lehmann et al.:
  # - destroy temporal order via shuffling (removes trends, autocorrelation)
  # - preserve spatial correlation (same permutation for all stations)
  # - recompute Ranom for each shuffled realization
  # - derive empirical CI from the ensemble
  #
  # ---- input:
  # df   : matrix/dataframe [stations x years], annual maxima (AMS) or similar
  # year : vector of observation years (length = ncol(df))
  # runs : number of Monte Carlo simulations (default = 1000)
  #
  # ---- output:
  # dataframe with:
  # - observed Ranom
  # - LOESS-smoothed trend
  # - 90% and 95% confidence intervals under iid null model
  #
    
  # length of time series
  index = seq(1, length(year), 1)
  
  # matrix to store simulated anomaly time series
  # rows = simulation runs
  # cols = time steps (sequence position, not calendar year)
  Ranom_multi <- matrix(NA, nrow = runs, ncol = length(year))
  
  # ---- Monte Carlo simulation loop
  for (i in 1:runs) {
    # randomly permute time indices (iid null hypothesis)
    # this destroys temporal structure (trend, persistence)
    index_sample = sample(index)

    # apply the same permutation to all stations
    # -> preserves spatial correlation between stations
    df_shuffle = df[, index_sample]

    # define synthetic time axis (1, 2, ..., N)
    # IMPORTANT:
    # after shuffling, "year" loses physical meaning
    # record statistics depend only on sequence order
    year_idx <- seq_along(year)

    # compute record-breaking anomaly on shuffled data
    # this yields Ranom under iid assumption
    df_anom <- RBF(df_shuffle, year = year_idx, plot = FALSE)

    # store anomaly time series for this simulation
    # each column corresponds to a sequence position
    Ranom_multi[i, ] <- df_anom$Ranom
  }
  
  # ---- Compute confidence intervals
  # For each time step (column), derive empirical quantiles

  # 90% CI (5th–95th percentile)
  ci_lower_90 <- apply(Ranom_multi, 2, quantile, probs = 0.05, na.rm = TRUE)
  ci_upper_90 <- apply(Ranom_multi, 2, quantile, probs = 0.95, na.rm = TRUE)
  
  # 95% CI (2.5th–97.5th percentile)
  ci_lower_95 <- apply(Ranom_multi, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_upper_95 <- apply(Ranom_multi, 2, quantile, probs = 0.975, na.rm = TRUE)

  # ---- Compute observed anomaly
  # This is the real-world Ranom time series (with actual years)
  df_anom = RBF(df, year, plot = FALSE)

  # ---- Compute smoothed trend (LOESS) 
  # NOTE:
  # Lehmann et al. use SSA; here LOESS is an approximation
  df_loess = data.frame(value = df_anom$Ranom, time = year)
  original_trend <- calculate_loess_trend(df_loess)
  
  # ---- Combine results 
  # Map simulation-based CI (sequence position) to actual years
  # (valid because time step index corresponds to progression in series)
  result = data.frame(
    time = year,                     # actual calendar year
    value = df_anom$Ranom,           # observed anomaly
    loess_trend = original_trend,    # smoothed trend
    ci_lower_90 = ci_lower_90,
    ci_upper_90 = ci_upper_90,
    ci_lower_95 = ci_lower_95,
    ci_upper_95 = ci_upper_95
  )
  return(result)
}


RB_Ranom_CI_index <- function(df, year, runs = 1000){
  # ---- description:
  # compute the confidence interval of record-breaking anomaly
  # through multiple shuffling and simulations
  # ---- input:
  # 
  # ---- output:
  # 
  # 
    
  index = seq(1, length(year), 1)
  
  Ranom_multi <- matrix(NA, nrow = runs, ncol = length(year))

  for (i in 1:runs) {
    index_sample = sample(index)
    df_shuffle = df[, index_sample]
    # bring estimated Ranom back to the original position
    df_anom = RBF(df_shuffle, year[index_sample], plot = FALSE)
    df_anom = df_anom[order(df_anom$year, decreasing = FALSE), ]
    
    Ranom_multi[i, ] <- df_anom$Ranom
  }
  
  ci_lower_90 <- apply(Ranom_multi, 2, quantile, probs = 0.05, na.rm = TRUE)
  ci_upper_90 <- apply(Ranom_multi, 2, quantile, probs = 0.95, na.rm = TRUE)
  
  ci_lower_95 <- apply(Ranom_multi, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_upper_95 <- apply(Ranom_multi, 2, quantile, probs = 0.975, na.rm = TRUE)

  df_anom = RBF(df, year, plot = FALSE)
  df_loess = data.frame(value = df_anom$Ranom, time = year)
  original_trend <- calculate_loess_trend(df_loess)
  
  result = data.frame(
    time = year,
    value = df_anom$Ranom,
    loess_trend = original_trend,
    ci_lower_90 = ci_lower_90,
    ci_upper_90 = ci_upper_90,
    ci_lower_95 = ci_lower_95,
    ci_upper_95 = ci_upper_95
  )
  return(result)
}




# ------ LOESS 

calculate_loess_trend <- function(data, span = 0.8) {
  loess_model <- loess(value ~ time, data = data, span = span)
  predict(loess_model, data.frame(time = data$time))
}

LOESS_Trend_CI <- function(df, year, runs = 1000){
  index = seq(1, length(year), 1)
  
  loess_trends <- matrix(NA, nrow = runs, ncol = length(year))
  
  for (i in 1:runs) {
    index_sample = sample(index)
    df_shuffle = df[, index_sample]
    df_anom = RBF(df_shuffle, year, plot = FALSE)
    
    df_loess = data.frame(value = df_anom$Ranom, time = year)
    
    loess_trends[i, ] <- calculate_loess_trend(df_loess)
  }
  
  ci_lower <- apply(loess_trends, 2, quantile, probs = 0.05, na.rm = TRUE)
  ci_upper <- apply(loess_trends, 2, quantile, probs = 0.95, na.rm = TRUE)
  
  df_anom = RBF(df, year, plot = FALSE)
  df_loess = data.frame(value = df_anom$Ranom, time = year)
  original_trend <- calculate_loess_trend(df_loess)
  
  plot_data = data.frame(
    time = year,
    value = df_anom$Ranom,
    loess_trend = original_trend,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  return(plot_data)
}


