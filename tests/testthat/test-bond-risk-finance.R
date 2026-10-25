.manual_bond_risk <- function(face, c, n, k, j, R = face) {
  N <- as.integer(round(n * k))
  coupon <- face * c / k
  period <- c(seq_len(N), N)
  cf <- c(rep(coupon, N), R)
  pv <- cf / (1 + j)^period
  price <- sum(pv)

  D_periods <- sum(period * pv) / price
  D_years <- D_periods / k
  C_periods <- sum(cf * period * (period + 1) / (1 + j)^(period + 2)) / price
  C_years <- C_periods / k^2

  list(
    price = price,
    D_periods = D_periods,
    D_years = D_years,
    C_periods = C_periods,
    C_years = C_years
  )
}

test_that("zero-coupon bond has Macaulay duration equal to maturity", {
  out <- bond_duration(
    face = 1000,
    c = 0,
    n = 5,
    k = 1,
    y_effective_per_period = 0.06
  )

  expect_equal(out$macaulay_duration_years, 5, tolerance = 1e-12)
  expect_equal(out$modified_duration_periods_j, 5 / 1.06, tolerance = 1e-12)
  expect_equal(out$modified_duration_years_i, 5 / 1.06, tolerance = 1e-12)
})

test_that("bond_duration agrees with direct discounted-cash-flow calculation", {
  expected <- .manual_bond_risk(
    face = 1000,
    c = 0.07,
    n = 8,
    k = 2,
    j = 0.025
  )

  out <- bond_duration(
    face = 1000,
    c = 0.07,
    n = 8,
    k = 2,
    y_effective_per_period = 0.025
  )

  expect_equal(out$price, expected$price, tolerance = 1e-10)
  expect_equal(out$macaulay_duration_periods, expected$D_periods, tolerance = 1e-10)
  expect_equal(out$macaulay_duration_years, expected$D_years, tolerance = 1e-10)
})

test_that("coupon bond Macaulay duration is positive and below maturity", {
  out <- bond_duration(
    face = 1000,
    c = 0.06,
    n = 10,
    k = 2,
    y_effective_per_period = 0.03
  )

  expect_gt(out$macaulay_duration_years, 0)
  expect_lt(out$macaulay_duration_years, 10)
})

test_that("zero-coupon convexity follows the discrete actuarial formula", {
  y <- 0.05
  n <- 6
  out <- bond_convexity(
    face = 1000,
    c = 0,
    n = n,
    k = 1,
    y_effective_per_period = y
  )

  expected <- n * (n + 1) / (1 + y)^2
  expect_equal(out$discrete_convexity_years, expected, tolerance = 1e-12)
})

test_that("bond_convexity agrees with direct discounted-cash-flow calculation", {
  expected <- .manual_bond_risk(
    face = 1000,
    c = 0.07,
    n = 8,
    k = 2,
    j = 0.025
  )

  out <- bond_convexity(
    face = 1000,
    c = 0.07,
    n = 8,
    k = 2,
    y_effective_per_period = 0.025
  )

  expect_equal(out$price, expected$price, tolerance = 1e-10)
  expect_equal(out$discrete_convexity_periods, expected$C_periods, tolerance = 1e-10)
  expect_equal(out$discrete_convexity_years, expected$C_years, tolerance = 1e-10)
})

test_that("bond duration and convexity use the same price basis", {
  dur <- bond_duration(
    face = 1000, c = 0.055, n = 7, k = 2,
    y_effective_per_period = 0.028
  )
  conv <- bond_convexity(
    face = 1000, c = 0.055, n = 7, k = 2,
    y_effective_per_period = 0.028
  )

  expect_equal(dur$price, conv$price, tolerance = 1e-12)
})

test_that("zero-term bond risk measures are zero", {
  dur <- bond_duration(
    face = 1000, c = 0.05, n = 0, k = 2,
    y_effective_per_period = 0.03, R = 1010
  )
  conv <- bond_convexity(
    face = 1000, c = 0.05, n = 0, k = 2,
    y_effective_per_period = 0.03, R = 1010
  )

  expect_equal(dur$price, 1010)
  expect_equal(dur$macaulay_duration_years, 0)
  expect_equal(conv$price, 1010)
  expect_equal(conv$discrete_convexity_years, 0)
})
