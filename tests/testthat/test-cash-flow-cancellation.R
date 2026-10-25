test_that("equal financial factors cancel before overflowing individual values", {
  for (order in list(1:3, c(3, 1, 2), c(2, 3, 1))) {
    cf <- c(2^1000, -2^1000, 7)[order]
    times <- c(1, 1, 0)[order]
    expect_equal(pv_flow(cf, i = -100, i_type = "force", t = times), 7)
    # At horizon 1, only the two opposite time-zero amounts accumulate.
    expect_equal(fv_flow(cf, i = 100, i_type = "force", t = 1 - times), 7)
  }
  expect_equal(pv_flow(c(1, -1), i = -800, i_type = "force", t = c(1, 1)), 0)
  expect_equal(fv_flow(c(1, -1, 0), i = 800, i_type = "force", t = c(0, 0, 1)), 0)
})

test_that("unequal times can share an exactly equal log factor", {
  # -delta_k * t_k = 100 for both of the opposite large amounts.
  expect_equal(pv_flow(c(2^1000, -2^1000, 9),
    i = c(-100, -50, 0), i_type = "force", t = c(1, 2, 0)), 9)
})

test_that("common-factor net amounts preserve representable residuals", {
  # (2^1000 - (2^1000 - 2^948)) * 2^50 = 2^998.
  # All amounts and the independent reference are exact binary powers.
  cf <- c(2^1000, -(2^1000 - 2^948))
  expect_equal(pv_flow(cf, i = -50 * log(2), i_type = "force", t = c(1, 1)) /
    2^998, 1, tolerance = 1e-12)
  expect_equal(fv_flow(c(cf, 0), i = 50 * log(2), i_type = "force", t = c(0, 0, 1)) /
    2^998, 1, tolerance = 1e-12)
})

test_that("grouping does not overflow a representable discounted total", {
  # Raw amounts sum to 2^1024 (outside double); discounted sum is 2^1022.
  expect_equal(pv_flow(c(2^1023, 2^1023), i = 3, t = c(1, 1)) / 2^1022, 1)
  expect_equal(pv_flow(c(1, 2, -1), i = 0.1, t = c(0, 1, 1)), 1 + 1 / 1.1)
  expect_identical(pv_flow(2^1000, i = -100, i_type = "force", t = 1), Inf)
  # Distinct extreme factors must not be falsely grouped by rounded labels.
  expect_true(is.nan(pv_flow(c(1, -1), i = c(-800, -801),
    i_type = "force", t = c(1, 1))))
})

test_that("plot present value shares common-factor cancellation", {
  p <- plot_cash_flow(C = c(2^1000, -2^1000, 7), t = c(1, 1, 0),
    i = -100, i_type = "force")
  expect_s3_class(p, "ggplot")
  expect_match(as.character(p$labels$subtitle), "7")
})
