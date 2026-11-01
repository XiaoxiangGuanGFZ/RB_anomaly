test_that("public compatibility aliases retain the exact implementations", {
  aliases <- c(
    Func_RB_matrix = "rb_record_matrix",
    Func_Ranom = "rb_record_indicator",
    Func_Ranom_1N = "rb_record_probability",
    Func_RB_region = "rb_anomaly_region",
    Func_RB_region_list = "rb_anomaly_region_list",
    ssa_trend = "rb_ssa_trend",
    ssa_trend_matrix = "rb_ssa_trend_matrix",
    calculate_loess_trend = "rb_loess_trend",
    senslope = "rb_sen_slope",
    lmslope = "rb_linear_slope",
    mk_test = "rb_mann_kendall",
    plot_Ranom = "rb_plot_anomaly",
    plot_Ranom_CI = "rb_plot_anomaly_ci",
    trend_formatter = "rb_format_trend",
    pvalue_formatter = "rb_format_pvalue",
    Func_RB_region_list_CI = "rb_legacy_raw_ci_list",
    RB_Ranom_CI = "rb_legacy_raw_ci",
    Func_RB_region_list_CI_fast = "rb_legacy_ssa_ci",
    LOESS_Trend_CI = "rb_legacy_loess_ci",
    RB_Ranom_CI_index = "rb_legacy_reindexed_ci"
  )
  exports <- getNamespaceExports("RBanomaly")
  for (old in names(aliases)) {
    expect_true(old %in% exports, info = old)
    expect_identical(getExportedValue("RBanomaly", old),
                     getFromNamespace(aliases[[old]], "RBanomaly"), info = old)
  }
  canonical <- unname(aliases[!grepl("^rb_legacy_", aliases)])
  expect_true(all(canonical %in% exports))
  expect_false(any(grepl("^rb_legacy_", exports)))
  expect_false(any(c("rb_validate_matrix", "rb_scan_records",
                     "rb_aggregate_records") %in% exports))
})
