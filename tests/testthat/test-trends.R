test_that("SSA reconstructs each row and preserves missing edges", {
  x <- rbind(c(NA, 1:20, NA), c(NA, 21:40, NA))
  fit <- rb_ssa_trend_matrix(x, L = 5)
  expect_equal(dim(fit), dim(x))
  expected <- Rssa::reconstruct(Rssa::ssa(1:20, L = 5, kind = "1d-ssa"),
                               groups = list(trend = 1))$trend
  expect_equal(fit[1, 2:21], as.numeric(expected))
  expect_true(all(is.na(fit[, c(1, 22)])))
  expect_equal(dim(rb_ssa_trend_matrix(x[1, , drop = FALSE], L = 5)), c(1L, 22L))
  expect_error(rb_ssa_trend(c(1, 2, NA, 4, 5), L = 2), "interior")
  expect_error(rb_ssa_trend(1:5, L = 5), "usable series length")
})

test_that("SSA Monte Carlo is reproducible and returns aligned intervals", {
  set.seed(12)
  x <- matrix(rexp(120), nrow = 4)
  a <- Func_RB_region_list_CI_fast(list(x), 1991:2020, runs = 4, L = 5)
  b <- Func_RB_region_list_CI_fast(list(x), 1991:2020, runs = 4, L = 5)
  expect_equal(a, b)
  expect_equal(a$time, 1991:2020)
  expect_true(all(is.na(a[1, -1])))
  expect_true(all(is.finite(as.matrix(a[-1, ]))))
  expect_true(all(a$ci_lower_95[-1] <= a$ci_lower_90[-1]))
  expect_true(all(a$ci_upper_95[-1] >= a$ci_upper_90[-1]))
  single <- Func_RB_region_list_CI_fast(list(x), 1991:2020, runs = 1, L = 5)
  expect_equal(nrow(single), 30L)
})

test_that("legacy confidence intervals work for a single station", {
  set.seed(5)
  x <- matrix(rexp(30), nrow = 1)
  for (fun in list(RB_Ranom_CI, RB_Ranom_CI_index, LOESS_Trend_CI)) {
    if (identical(fun, RB_Ranom_CI_index)) {
      expect_warning(result <- fun(x, 1991:2020, runs = 3), "deprecated")
    } else result <- fun(x, 1991:2020, runs = 3)
    expect_equal(nrow(result), 30L)
    expect_equal(result$time, 1991:2020)
    expect_true(all(is.finite(result$loess_trend[-1])))
  }
  a <- Func_RB_region_list_CI(list(x), 1991:2020, runs = 3)
  b <- Func_RB_region_list_CI(list(x), 1991:2020, runs = 3)
  expect_equal(a, b)
  expect_equal(nrow(a), 30L)
})
