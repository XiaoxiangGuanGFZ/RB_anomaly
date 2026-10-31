test_that("example data load through the installed package interface", {
  e <- new.env(parent = emptyenv())
  data("mega_floods", package = "RBanomaly", envir = e)
  expect_equal(dim(e$df_AMS_All), c(218L, 8024L))
  expect_equal(dim(e$df_gauge), c(8023L, 8L))
})

test_that("plot helpers create renderable ggplot objects", {
  p <- plot_Ranom(c(0, 20, -10), 2001:2003, p_print = FALSE)
  expect_s3_class(p, "ggplot")
  expect_length(ggplot2::ggplot_build(p)$data, 2L)
  df <- data.frame(time = 1:3, value = c(0, 20, -10), loess_trend = 0,
                   ci_lower_90 = -20, ci_upper_90 = 20,
                   ci_lower_95 = -30, ci_upper_95 = 30)
  p <- plot_Ranom_CI(df, p_print = FALSE)
  expect_length(ggplot2::ggplot_build(p)$data, 5L)
})

test_that("trend and table helpers resolve their namespace dependencies", {
  expect_equal(unname(lmslope(1:20)), 1)
  expect_length(senslope(1:20), 2L)
  expect_length(mk_test(1:20), 2L)
  expect_length(trend_formatter(c(-1, 1)), 2L)
  expect_length(pvalue_formatter(c(0.01, 0.5)), 2L)
})
