
library(tidyverse)
library(trend)
library(ggplot2)


calculate_loess_trend <- function(data, span = 0.8) {
  loess_model = loess(value ~ time, data = data, span = span)
  predict(loess_model, data.frame(time = data$time))
}

# ----------- Record-breaking anomaly: matrix-based computation ----------

Func_RB_matrix <- function(df) {
  # ============================================================
  # Function: Func_RB_matrix
  # ============================================================
  # Description:
  # Fast matrix-based computation of:
  #   (1) observed record-breaking events (Robs)
  #   (2) expected record-breaking probability under iid climate (R1n)
  #
  # Compared to the original implementation:
  # - avoids repeatedly scanning past values
  # - updates running maximum incrementally
  # - computes all stations simultaneously
  #
  # This reduces complexity from roughly O(n^2) per station
  # to approximately O(n).
  #
  # ------------------------------------------------------------
  # Input:
  # df : matrix/data.frame
  #      rows    = stations/gauges
  #      columns = time (years)
  #
  # ------------------------------------------------------------
  # Output:
  # A list with:
  #   Robs : matrix of observed records
  #           1 = record-breaking event
  #           0 = not a record
  #          NA = invalid / insufficient history
  #
  #   R1n  : matrix of expected record probability under iid:
  #           1 / N
  #           where N = number of valid observations so far
  #
  # ------------------------------------------------------------
  # Notes:
  # - Missing values are allowed.
  # - Spatial aggregation is NOT performed here.
  # - This function only computes station-level quantities.
  # ============================================================
  
  # ---- convert input to numeric matrix
  x <- as.matrix(df)
  x <- apply(x, 2, as.numeric)
  x <- matrix(x, nrow = nrow(df))

  # dimensions
  p <- nrow(x)   # number of stations
  n <- ncol(x)   # number of time steps

  # ---- storage for results
  Robs <- matrix(NA_real_, p, n)
  R1n  <- matrix(NA_real_, p, n)

  # ----------------------------------------------------------
  # running statistics
  # ----------------------------------------------------------
  # past_max:
  # current maximum value observed so far at each station
  #
  # past_n:
  # number of valid observations encountered so far
  #
  # Both vectors are updated sequentially in time.
  # ----------------------------------------------------------

  past_max <- rep(-Inf, p)
  past_n   <- rep(0, p)

  # ==========================================================
  # iterate over time
  # ==========================================================
  for (j in 1:n) {

    # stations with valid value at current time
    valid_now <- !is.na(x[, j])

    # stations with at least one previous valid observation
    has_past <- past_n > 0

    # stations eligible for RB evaluation
    use <- valid_now & has_past

    # --------------------------------------------------------
    # observed record-breaking event
    #
    # x > previous running maximum
    # --------------------------------------------------------
    Robs[use, j] <- as.numeric(
      x[use, j] > past_max[use]
    )

    # --------------------------------------------------------
    # expected record probability under iid assumption
    #
    # probability = 1 / N
    # where N = number of valid observations including current
    # --------------------------------------------------------
    R1n[use, j] <- 1 / (past_n[use] + 1)

    # --------------------------------------------------------
    # update running history AFTER evaluating current year
    # --------------------------------------------------------
    update <- valid_now

    # update running maximum
    past_max[update] <- pmax(
      past_max[update],
      x[update, j]
    )

    # update number of valid observations
    past_n[update] <- past_n[update] + 1
  }

  return(list(
    Robs = Robs,
    R1n  = R1n
  ))

  # The pmax() function in R obtains the parallel maxima of 
  # two or more given vectors. In simpler words, 
  # the function compares the values of the two vectors in pairs and 
  # returns a final vector that contains the maximum values from each pair.

}


Func_RB_region <- function(df, year, plot = FALSE) {
  # ============================================================
  # Function: Func_RB_region
  # ============================================================
  # Description:
  # Compute regional record-breaking anomaly (Ranom)
  # across multiple stations/gauges.
  #
  # The function:
  #   1. computes station-level record-breaking statistics
  #   2. aggregates them spatially across stations
  #   3. computes regional anomaly:
  #
  #      Ranom = (Robs - R1n) / R1n * 100
  #
  # ------------------------------------------------------------
  # Input:
  # df   : matrix/data.frame
  #        rows    = stations
  #        columns = years
  #
  # year : vector of years
  #
  # plot : logical
  #        TRUE  -> plot anomaly series
  #        FALSE -> no plot
  #
  # ------------------------------------------------------------
  # Output:
  # data.frame with:
  #   year
  #   R1n   : expected records under iid climate
  #   Robs  : observed records
  #   Ranom : record-breaking anomaly (%)
  # ============================================================
  # ---- consistency check
  if (length(year) != ncol(df)) {
    stop(
      "Stopped! Length of year doesn't match number of columns in df."
    )
  }

  # ----------------------------------------------------------
  # compute station-level record statistics
  # ----------------------------------------------------------
  rb <- Func_RB_matrix(df)

  # ----------------------------------------------------------
  # aggregate spatially across stations
  # ----------------------------------------------------------
  out <- data.frame(
    year = year,

    # expected record counts
    R1n = colSums(rb$R1n, na.rm = TRUE),

    # observed record counts
    Robs = colSums(rb$Robs, na.rm = TRUE)
  )

  # ----------------------------------------------------------
  # compute RB anomaly (%)
  # ----------------------------------------------------------
  out$Ranom <- ifelse(
    out$R1n > 0,
    round(
      (out$Robs - out$R1n) / out$R1n * 100,
      1
    ),
    NA
  )

  # ----------------------------------------------------------
  # optional visualization
  # ----------------------------------------------------------
  if (plot) {

    p <- ggplot(out) +
      geom_col(aes(x = year, y = Ranom)) +
      labs(
        x = "time",
        y = "record-breaking anomaly (%)"
      ) +
      theme_minimal()

    print(p)
  }

  return(out)
}



Func_RB_region_list <- function(
    list_df,
    year,
    plot = FALSE
) {
  # ============================================================
  # Function: Func_RB_region_list
  # ============================================================
  # Description:
  # Compute record-breaking anomaly for a LIST of datasets
  # and aggregate them into one combined anomaly series.
  #
  # Typical application:
  # - each list element represents one month
  # - compute monthly RB anomalies
  # - aggregate to seasonal or annual scale
  #
  # ------------------------------------------------------------
  # Input:
  # list_df : list of matrices/data.frames
  #
  #           each element:
  #             rows    = stations
  #             columns = years
  #
  # year    : vector of years
  #
  # plot    : logical
  #
  # ------------------------------------------------------------
  # Output:
  # A list containing:
  #
  #   [[1]], [[2]], ...
  #       regional RB anomaly for each element
  #
  #   $out_Ranom
  #       aggregated RB anomaly across all elements
  #
  # ------------------------------------------------------------
  # Notes:
  # - aggregation is performed using SUMS (not means)
  # - this follows the Lehmann et al. framework
  # ============================================================

  # number of list elements
  n_list <- length(list_df)

  # ----------------------------------------------------------
  # ensure all elements have identical dimensions
  # ----------------------------------------------------------
  dims <- lapply(list_df, dim)

  if (!all(
    vapply(dims, identical, logical(1), dims[[1]])
  )) {
    stop(
      "RB anomaly calculation stopped. Elements have different shapes."
    )
  }

  # number of time steps
  time_len <- ncol(list_df[[1]])

  # ----------------------------------------------------------
  # storage for aggregated monthly/seasonal results
  # ----------------------------------------------------------
  Robs <- matrix(
    NA_real_,
    nrow = n_list,
    ncol = time_len
  )

  R1n <- matrix(
    NA_real_,
    nrow = n_list,
    ncol = time_len
  )

  # storage for outputs
  output_list <- vector("list", n_list)

  # ==========================================================
  # iterate over list elements
  # ==========================================================
  for (i in seq_len(n_list)) {

    # regional RB anomaly for this element
    out_i <- Func_RB_region(
      list_df[[i]],
      year,
      plot = FALSE
    )

    # store full output
    output_list[[i]] <- out_i

    # extract Robs and R1n for aggregation
    Robs[i, ] <- out_i$Robs
    R1n[i, ]  <- out_i$R1n
  }

  # ----------------------------------------------------------
  # aggregate across all list elements
  # ----------------------------------------------------------
  out_Ranom <- data.frame(
    time = year,

    # aggregate observed records
    Robs = colSums(Robs, na.rm = TRUE),

    # aggregate expected records
    R1n = colSums(R1n, na.rm = TRUE)
  )

  # ----------------------------------------------------------
  # compute final anomaly (%)
  # ----------------------------------------------------------
  out_Ranom$Ranom <- ifelse(
    out_Ranom$R1n > 0,
    round(
      (out_Ranom$Robs - out_Ranom$R1n) /
        out_Ranom$R1n * 100,
      1
    ),
    NA
  )

  # attach aggregated result
  output_list$out_Ranom <- out_Ranom

  # ----------------------------------------------------------
  # optional visualization
  # ----------------------------------------------------------
  if (plot) {

    p <- ggplot(out_Ranom) +
      geom_col(aes(x = time, y = Ranom)) +
      labs(
        x = "time",
        y = "record-breaking anomaly (%)"
      ) +
      theme_minimal()

    print(p)
  }

  return(output_list)
}



# -------------- RB anomaly: CI with multiple simulations ---------------------

Func_prepare_matrix_list <- function(list_df) {
  
    # ============================================================
    # Function: Func_prepare_matrix_list
    # Helper: convert list elements to numeric matrices once
    # ============================================================
    # Description:
    # Convert all elements in a list of data frames/matrices
    # into numeric matrices.
    #
    # This preprocessing step is important for performance:
    # - avoids repeated coercion during Monte Carlo runs
    # - ensures fast matrix operations
    # - guarantees numeric storage mode
    #
    # ------------------------------------------------------------
    # Input:
    # list_df : list of data.frames or matrices
    #
    # ------------------------------------------------------------
    # Output:
    # A list of numeric matrices with the same dimensions
    # as the original inputs.
    # ============================================================
  lapply(list_df, function(df) {

    # convert to matrix
    x <- as.matrix(df)

    # ensure numeric type
    x <- apply(x, 2, as.numeric)

    # restore matrix shape
    matrix(x, nrow = nrow(df))
  })
}



Func_RB_matrix_sums <- function(x, perm = NULL) {
  
    # ============================================================
    # Function: Func_RB_matrix_sums
    # Helper: compute regional Robs and R1n sums directly
    # ============================================================
    # Description:
    # Fast computation of spatially aggregated:
    #
    #   (1) observed record-breaking events (Robs)
    #   (2) expected iid record counts (R1n)
    #
    # for ALL stations simultaneously.
    #
    # Unlike earlier implementations:
    # - does NOT store station-level matrices
    # - computes only aggregated sums
    # - designed specifically for Monte Carlo simulations
    #
    # This substantially reduces memory allocation and runtime.
    #
    # ------------------------------------------------------------
    # Input:
    # x    : numeric matrix
    #        rows    = stations
    #        columns = years
    #
    # perm : optional permutation vector
    #
    #        NULL:
    #          use original temporal order
    #
    #        integer vector:
    #          shuffled temporal order used for iid null model
    #
    # ------------------------------------------------------------
    # Output:
    # A list with:
    #
    #   Robs : vector of aggregated observed records
    #   R1n  : vector of aggregated expected records
    #
    # Both vectors have length = number of years.
    #
    # ------------------------------------------------------------
    # Notes:
    # - the same temporal permutation should be applied
    #   to all variables/months to preserve spatial correlation
    #
    # - missing values are allowed
    #
    # - record statistics depend only on sequence order,
    #   not on actual calendar years
    # ============================================================

  # number of stations
  p <- nrow(x)

  # number of years
  n <- ncol(x)

  # ----------------------------------------------------------
  # if no permutation supplied:
  # use ORIGINAL order
  # ----------------------------------------------------------
  if (is.null(perm)) {
    perm <- seq_len(n)
  }

  # ----------------------------------------------------------
  # storage for spatially aggregated results
  # ----------------------------------------------------------
  Robs_sum <- rep(0, n)
  R1n_sum  <- rep(0, n)

  # ----------------------------------------------------------
  # running station-wise history
  #
  # past_max:
  #   current maximum value encountered so far
  #
  # past_n:
  #   number of valid observations encountered so far
  # ----------------------------------------------------------
  past_max <- rep(-Inf, p)
  past_n   <- rep(0, p)

  # ==========================================================
  # iterate over temporal sequence
  # ==========================================================
  for (j in seq_len(n)) {

    # actual column index after permutation
    col_j <- perm[j]

    # values at current time step
    x_now <- x[, col_j]

    # stations with valid current observation
    valid_now <- !is.na(x_now)

    # stations with at least one previous valid observation
    has_past <- past_n > 0

    # stations eligible for RB computation
    use <- valid_now & has_past

    # --------------------------------------------------------
    # compute observed and expected records
    # only for eligible stations
    # --------------------------------------------------------
    if (any(use)) {

      # observed record-breaking events
      Robs_sum[j] <- sum(
        x_now[use] > past_max[use]
      )

      # expected iid records
      R1n_sum[j] <- sum(
        1 / (past_n[use] + 1)
      )
    }

    # --------------------------------------------------------
    # update running history AFTER evaluation
    # --------------------------------------------------------
    update <- valid_now

    # update running maxima
    past_max[update] <- pmax(
      past_max[update],
      x_now[update]
    )

    # update number of valid observations
    past_n[update] <- past_n[update] + 1
  }

  return(list(
    Robs = Robs_sum,
    R1n  = R1n_sum
  ))
}




Func_RB_region_list_Ranom_only <- function(list_mat, perm = NULL) {
  
    # ============================================================
    # Function: Func_RB_region_list_Ranom_only
    # Helper: compute aggregated Ranom for a list of matrices
    # ============================================================
    # Description:
    # Compute aggregated regional record-breaking anomaly
    # for a LIST of matrices.
    #
    # This is a lightweight implementation specifically
    # optimized for Monte Carlo simulations.
    #
    # Compared to the full workflow:
    # - computes ONLY the quantities needed for CI estimation
    # - avoids creation of intermediate data frames/lists
    # - minimizes memory allocation
    #
    # ------------------------------------------------------------
    # Typical application:
    #
    # list_mat[[1]] = January AMS matrix
    # list_mat[[2]] = February AMS matrix
    # ...
    #
    # Aggregate all months to obtain annual/seasonal anomaly.
    #
    # ------------------------------------------------------------
    # Input:
    # list_mat : list of numeric matrices
    #
    # perm     : optional permutation vector
    #
    #            NULL:
    #              original temporal order
    #
    #            integer vector:
    #              shuffled temporal order
    #
    # ------------------------------------------------------------
    # Output:
    # data.frame with:
    #
    #   Robs  : aggregated observed records
    #   R1n   : aggregated expected iid records
    #   Ranom : record-breaking anomaly (%)
    #
    # ------------------------------------------------------------
    # Notes:
    # - aggregation uses SUMS (not means)
    # - follows Lehmann et al. methodology
    # ============================================================

  # number of years
  n <- ncol(list_mat[[1]])

  # ----------------------------------------------------------
  # storage for total aggregated records
  # ----------------------------------------------------------
  Robs_total <- rep(0, n)
  R1n_total  <- rep(0, n)

  # ==========================================================
  # iterate over list elements
  # ==========================================================
  for (i in seq_along(list_mat)) {

    # compute aggregated records for this element
    out_i <- Func_RB_matrix_sums(
      list_mat[[i]],
      perm = perm
    )

    # accumulate totals
    Robs_total <- Robs_total + out_i$Robs
    R1n_total  <- R1n_total  + out_i$R1n
  }

  # ----------------------------------------------------------
  # compute record-breaking anomaly (%)
  # ----------------------------------------------------------
  Ranom <- ifelse(
    R1n_total > 0,
    (Robs_total - R1n_total) /
      R1n_total * 100,
    NA
  )

  return(data.frame(
    Robs  = Robs_total,
    R1n   = R1n_total,
    Ranom = Ranom
  ))
}



Func_RB_region_list_CI_fast <- function(
    list_df,
    year,
    seed = 1,
    runs = 1000
) {
  
    # ============================================================
    # Function: Func_RB_region_list_CI_fast
    # Fast Monte Carlo CI for record-breaking anomaly
    # ============================================================
    # Description:
    # Fast Monte Carlo estimation of confidence intervals (CI)
    # for regional record-breaking anomaly (Ranom).
    #
    # The method follows the iid null-model framework:
    #
    #   1. randomly shuffle temporal order
    #   2. preserve spatial correlation
    #   3. recompute record statistics
    #   4. derive empirical confidence intervals
    #
    # ------------------------------------------------------------
    # Performance improvements:
    #
    # - matrices are prepared once
    # - no repeated data-frame creation
    # - no repeated station-level outputs
    # - only aggregated Ranom is computed
    #
    # This version is substantially faster than the original
    # implementation and suitable for large Monte Carlo ensembles.
    #
    # ------------------------------------------------------------
    # Input:
    # list_df : list of matrices/data.frames
    #           rows    = stations
    #           columns = years
    #
    # year    : vector of observation years
    #
    # seed    : random seed
    #
    # runs    : number of Monte Carlo simulations
    #
    # ------------------------------------------------------------
    # Output:
    # data.frame with:
    #
    #   time
    #   value           : observed Ranom
    #   loess_trend     : smoothed trend
    #
    #   ci_lower_90
    #   ci_upper_90
    #
    #   ci_lower_95
    #   ci_upper_95
    #
    # ------------------------------------------------------------
    # Notes:
    # - confidence intervals are derived empirically
    # - temporal order is destroyed during shuffling
    # - actual years lose physical meaning in simulations
    # - sequence position alone determines record probability
    # ============================================================
  
  # reproducibility
  set.seed(seed)

  # ----------------------------------------------------------
  # convert all inputs to numeric matrices ONCE
  # ----------------------------------------------------------
  list_mat <- Func_prepare_matrix_list(list_df)

  # number of years
  n <- length(year)

  # ----------------------------------------------------------
  # consistency checks
  # ----------------------------------------------------------
  if (ncol(list_mat[[1]]) != n) {
    stop(
      "Length of year does not match number of columns."
    )
  }

  # all matrices must have identical dimensions
  dims <- lapply(list_mat, dim)

  if (!all(
    vapply(dims, identical, logical(1), dims[[1]])
  )) {
    stop(
      "Elements in list_df have different shapes."
    )
  }

  # ----------------------------------------------------------
  # storage for Monte Carlo anomaly series
  #
  # rows = simulation runs
  # cols = sequence positions
  # ----------------------------------------------------------
  Ranom_multi <- matrix(
    NA_real_,
    nrow = runs,
    ncol = n
  )

  # ==========================================================
  # Monte Carlo simulation loop
  # ==========================================================
  for (r in seq_len(runs)) {

    # --------------------------------------------------------
    # generate random temporal permutation
    #
    # IMPORTANT:
    # the SAME permutation is applied to all matrices
    # to preserve spatial and cross-variable correlation
    # --------------------------------------------------------
    perm <- sample.int(n)

    # compute shuffled anomaly
    sim_out <- Func_RB_region_list_Ranom_only(
      list_mat,
      perm = perm
    )

    # store anomaly series
    Ranom_multi[r, ] <- sim_out$Ranom
  }

  # ----------------------------------------------------------
  # empirical confidence intervals
  # ----------------------------------------------------------

  # 90% confidence interval
  ci_lower_90 <- apply(
    Ranom_multi,
    2,
    quantile,
    probs = 0.05,
    na.rm = TRUE
  )

  ci_upper_90 <- apply(
    Ranom_multi,
    2,
    quantile,
    probs = 0.95,
    na.rm = TRUE
  )

  # 95% confidence interval
  ci_lower_95 <- apply(
    Ranom_multi,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  )

  ci_upper_95 <- apply(
    Ranom_multi,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # compute observed anomaly (original temporal order)
  # ----------------------------------------------------------
  obs_out <- Func_RB_region_list_Ranom_only(
    list_mat,
    perm = NULL
  )

  # ----------------------------------------------------------
  # smooth observed anomaly using LOESS
  #
  # NOTE:
  # Lehmann et al. used SSA
  # here LOESS acts as approximation
  # ----------------------------------------------------------
  df_loess <- data.frame(
    value = obs_out$Ranom,
    time = year
  )

  original_trend <- calculate_loess_trend(df_loess)

  # ----------------------------------------------------------
  # final output
  # ----------------------------------------------------------
  return(data.frame(
    time = year,
    value = obs_out$Ranom,
    loess_trend = original_trend,

    ci_lower_90 = ci_lower_90,
    ci_upper_90 = ci_upper_90,

    ci_lower_95 = ci_lower_95,
    ci_upper_95 = ci_upper_95
  ))
}

