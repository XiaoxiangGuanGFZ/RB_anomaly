test_that("enumerated permutations recover the exact continuous iid record law", {
  permutations <- function(x) {
    if (length(x) == 1L) return(matrix(x, nrow = 1))
    do.call(rbind, lapply(seq_along(x), function(i) {
      cbind(x[i], permutations(x[-i]))
    }))
  }
  records <- rb_record_matrix(permutations(1:5))
  expect_equal(colMeans(records$Robs)[-1], 1 / (2:5))
  expect_equal(mean(rowSums(records$Robs, na.rm = TRUE)), sum(1 / (2:5)))
})

test_that("edge cases retain dimensions, eligibility, and exact count pooling", {
  x <- rbind(up = c(NA, 1, 2, 3, NA), down = c(NA, 3, 2, 1, NA),
             tied = c(NA, 2, 2, 2, NA), absent = rep(NA_real_, 5),
             single = c(NA, NA, 2, NA, NA))
  colnames(x) <- as.character(2001:2005)
  detail <- rb_record_matrix(x)
  out <- rb_anomaly_region(x, 2001:2005)
  expect_identical(dimnames(detail$Robs), dimnames(x))
  expect_equal(out$Robs, c(0, 0, 1, 1, 0))
  expect_equal(out$R1n, c(0, 0, 1.5, 1, 0))
  expect_equal(out$n_eligible, c(0, 0, 3, 3, 0))
  expect_equal(out$Ranom, c(NA, NA, -100/3, 0, NA))
  one <- rb_record_matrix(matrix(5, 1, 1))
  expect_equal(dim(one$Robs), c(1L, 1L))
  expect_true(is.na(one$Robs[1, 1]))
  expect_error(rb_record_indicator(c(NA, 1)), "two valid")
  pooled <- rb_anomaly_region_list(list(matrix(c(1, 2, 3), 1),
                                   matrix(c(NA, 2, 1), 1)), 1:3)$out_Ranom
  expect_equal(pooled$Ranom[3], 100 * (1 / (1/3 + 1/2) - 1))
  expect_equal(pooled$n_eligible, c(0, 1, 2))
})

test_that("input validation rejects ambiguous data before computation", {
  x <- matrix(1:30, 1)
  expect_error(rb_record_matrix(matrix(character(), 0, 0)), "non-empty")
  expect_error(rb_record_matrix(data.frame(a = "1")), "numeric")
  expect_error(rb_record_matrix(matrix(c(1, Inf), 1)), "finite")
  expect_error(rb_anomaly_region(x, rep(1, 30)), "strictly increasing")
  expect_error(rb_anomaly_region_list(list(), 1:30), "non-empty")
  expect_error(rb_anomaly_region_list(list(x, t(x)), 1:30), "shapes")
  expect_error(rb_anomaly_ci(list(x), 1:30, runs = 0), "positive integer")
  expect_error(rb_anomaly_ci(list(x), 1:30, runs = 1.5), "positive integer")
  expect_error(rb_anomaly_ci(list(x), 1:30, seed = -1), "seed")
  expect_error(rb_anomaly_ci(list(x), c(1:29, 31)), "regularly spaced")
  expect_error(rb_anomaly_ci(list(x), 1:30, L = 30), "usable")
  expect_error(rb_anomaly_ci(list(x), 1:30, statistic = "loess", span = 0), "span")
  x[1, 10] <- NA
  expect_error(rb_anomaly_ci(list(x), 1:30, runs = 1), "interior")
})

test_that("null bands match an independent manual permutation calculation", {
  x <- rbind(seq_len(30), rev(seq_len(30)), sin(seq_len(30)))
  y <- x * 2 + 1
  year <- 1991:2020
  # Manual record calculation deliberately scans previous observations.
  manual <- function(matrices) {
    obs <- exp <- numeric(ncol(matrices[[1]]))
    for (m in matrices) for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
      past <- m[i, seq_len(j - 1L)]
      past <- past[!is.na(past)]
      if (!is.na(m[i, j]) && length(past)) {
        obs[j] <- obs[j] + as.numeric(m[i, j] > max(past))
        exp[j] <- exp[j] + 1 / (length(past) + 1)
      }
    }
    ifelse(exp > 0, 100 * (obs - exp) / exp, NA_real_)
  }
  set.seed(73)
  permutations <- replicate(8, sample.int(30), simplify = FALSE)
  raw <- t(vapply(permutations, function(p) manual(list(x[, p], y[, p])), numeric(30)))
  for (target in c("raw", "ssa", "loess")) {
    smooth <- switch(target,
      raw = function(z) z,
      ssa = function(z) rb_ssa_trend(z, L = 5),
      loess = function(z) rb_loess_trend(data.frame(time = year, value = z)))
    simulations <- t(apply(raw, 1, smooth))
    expected <- apply(simulations[, -1, drop = FALSE], 2, quantile,
                      probs = c(.05, .95, .025, .975), names = FALSE)
    fit <- rb_anomaly_ci(list(x, y), year, target, runs = 8, seed = 73, L = 5)
    bands <- as.matrix(fit$data[-1, c("ci_lower_90", "ci_upper_90",
                                      "ci_lower_95", "ci_upper_95")])
    expect_equal(unname(bands), unname(t(expected)))
    expect_equal(fit$data$statistic_value, smooth(manual(list(x, y))))
    expect_identical(fit$settings$statistic, target)
    expect_equal(fit$data$time, year)
    expect_s3_class(rb_plot_anomaly_ci(fit, p_print = FALSE), "ggplot")
  }
})

test_that("seeded runs are local and wrappers retain their targets", {
  x <- matrix(seq_len(30), 1)
  set.seed(24)
  before <- .Random.seed
  fit <- rb_anomaly_ci(list(x), 1:30, "raw", runs = 4, seed = 9)
  expect_identical(.Random.seed, before)
  expect_equal(fit, rb_anomaly_ci(list(x), 1:30, "raw", runs = 4, seed = 9))
  legacy <- Func_RB_region_list_CI(list(x), 1:30, seed = 9, runs = 4)
  expect_equal(legacy$ci_lower_95, fit$data$ci_lower_95)
  set.seed(9)
  expect_equal(RB_Ranom_CI(x, 1:30, runs = 4), legacy)
  ssa <- rb_anomaly_ci(list(x), 1:30, runs = 4, seed = 9, L = 5)
  old <- Func_RB_region_list_CI_fast(list(x), 1:30, seed = 9, runs = 4, L = 5)
  expect_equal(old$loess_trend, ssa$data$statistic_value)
  expect_equal(old$ci_lower_95, ssa$data$ci_lower_95)
  before <- .Random.seed
  expect_error(rb_anomaly_ci(list(x), 1:30, L = 30), "usable")
  expect_identical(.Random.seed, before)
  set.seed(9)
  loess <- LOESS_Trend_CI(x, 1:30, runs = 4)
  new <- rb_anomaly_ci(list(x), 1:30, "loess", runs = 4, seed = 9)
  expect_equal(loess$ci_lower, new$data$ci_lower_90)
  expect_length(ggplot2::ggplot_build(rb_plot_anomaly_ci(loess, FALSE))$data, 5)
  expect_warning(RB_Ranom_CI_index(x, 1:30, runs = 2), "deprecated")
})

test_that("raw null bands handle unavailable years and single observations", {
  fit <- rb_anomaly_ci(list(matrix(NA_real_, 2, 4)), 1:4, "raw", runs = 1)
  expect_true(all(is.na(fit$data$ci_lower_95)))
  expect_equal(fit$data$n_eligible, rep(0, 4))
  one <- rb_anomaly_ci(list(matrix(1, 1, 1)), 2000, "raw", runs = 1)
  expect_true(is.na(one$data$statistic_value))
  expect_equal(one$data$Robs, 0)
})
