test_that("record indicators handle ties, missing years, and delayed starts", {
  x <- rbind(c(1, 2, 2, NA, 4), c(NA, 3, 1, NA, 5))
  records <- rb_record_matrix(x)
  expect_equal(records$Robs,
               rbind(c(NA, 1, 0, NA, 1), c(NA, NA, 0, NA, 1)))
  expect_equal(records$R1n,
               rbind(c(NA, 1/2, 1/3, NA, 1/4), c(NA, NA, 1/2, NA, 1/3)))
  expect_equal(rb_record_indicator(x[1, ]), records$Robs[1, ])
  expect_equal(rb_record_probability(x[1, ]), records$R1n[1, ])
  regional <- rb_anomaly_region(x, 2001:2005)
  expect_equal(regional$Robs, c(0, 1, 0, 0, 2))
  expect_equal(regional$Ranom, c(NA, 100, -100, NA, 100 * (2 / (1/4 + 1/3) - 1)))
  expect_error(rb_anomaly_region(x, 1:2), "Length of year")
})

test_that("list aggregation sums counts and supports a single dataset", {
  x <- rbind(1:10, 10:1)
  one <- rb_anomaly_region_list(list(x), 1:10)$out_Ranom
  two <- rb_anomaly_region_list(list(x, x), 1:10)$out_Ranom
  expect_equal(two$Robs, 2 * one$Robs)
  expect_equal(two$R1n, 2 * one$R1n)
  expect_equal(two$Ranom, one$Ranom)
})

test_that("fast aggregation agrees with the public regional result", {
  x <- rbind(c(NA, 1, 2, 2, NA, 3), c(2, 1, 4, 3, NA, 5))
  fast <- RBanomaly:::rb_anomaly_from_matrices(list(x))
  full <- rb_anomaly_region_list(list(x), 1:6)$out_Ranom
  usable <- !is.na(full$Ranom)
  expect_equal(fast$Robs[usable], full$Robs[usable])
  expect_equal(fast$R1n[usable], full$R1n[usable])
  expect_equal(fast$Ranom, full$Ranom)
  expect_equal(fast$Robs, full$Robs)
  expect_equal(fast$R1n, full$R1n)
  expect_true(all(is.na(fast$Ranom[c(1, 5)])))
})
