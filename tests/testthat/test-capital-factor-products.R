test_that("capital and factor are combined before declaring numerical loss", {
  for (sgn in c(-1, 1)) {
    for (power in c(-1000, 1000)) {
      capital <- sgn * 2^power
      force <- -sign(power) * 1100 * log(2)
      expected <- sgn * 2^(power - sign(power) * 1100)
      # Binary powers provide an independent reference: 2^a * 2^b=2^(a+b).
      expect_equal(future_value(C = capital, i = force, i_type = "force", t = 1) /
        expected, 1, tolerance = 1e-12)
      expect_equal(present_value(C = capital, i = -force, i_type = "force", t = 1) /
        expected, 1, tolerance = 1e-12)
      expect_equal(pv_flow(cf = capital, i = -force, i_type = "force", t = 1) /
        expected, 1, tolerance = 1e-12)
      expect_equal(fv_flow(cf = c(capital, 0), i = force, i_type = "force", t = c(0, 1)) /
        expected, 1, tolerance = 1e-12)
    }
  }
})

test_that("known zero payments remain zero with unrepresentable factors", {
  expect_equal(future_value(C = 0, i = 800, i_type = "force", t = 1), 0)
  expect_equal(present_value(C = 0, i = -800, i_type = "force", t = 1), 0)
  expect_equal(pv_flow(cf = c(0, 0), i = -800, i_type = "force", t = c(1, 2)), 0)
  expect_equal(fv_flow(cf = c(0, 0), i = 800, i_type = "force", t = c(0, 1)), 0)
  expect_true(is.na(present_value(C = 0, i = NA_real_, t = 1)))
  expect_true(is.na(future_value(C = 0, i = 0, t = NA_real_)))
})

test_that("a subnormal intermediate factor does not degrade a normal value", {
  # Decimal arithmetic at 80 digits, on exact inputs: 2^1000 * exp(-740).
  expected <- 4.4882708347358573677923533961594328384e-21
  expect_equal(present_value(C = 2^1000, i = 740, i_type = "force", t = 1) /
    expected, 1, tolerance = 2e-12)
})
