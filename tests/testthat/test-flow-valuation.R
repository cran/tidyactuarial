test_that("pv_flow agrees with direct discounting under a constant rate", {
  cf <- c(-1000, 300, 400, 500)
  t <- c(0, 1, 2, 3)
  i <- 0.08

  expected <- sum(cf / (1 + i)^t)

  expect_equal(
    pv_flow(cf = cf, i = i, t = t),
    expected
  )
})


test_that("pv_flow discounts each cash flow with its associated spot rate", {
  cf <- c(100, 150, 200)
  t <- c(1, 2, 3)
  spot <- c(0.05, 0.055, 0.06)

  expected <- sum(cf / (1 + spot)^t)

  expect_equal(
    pv_flow(
      cf = cf,
      i = spot,
      i_type = "effective",
      t = t
    ),
    expected
  )
})


test_that("pv_flow accepts heterogeneous rate conventions", {
  cf <- c(100, 150, 200)
  t <- c(1, 2, 3)
  rates <- c(0.05, 0.12, log(1.07))
  types <- c("effective", "nominal_interest", "force")
  m <- c(1, 12, 1)

  i_eff <- standardize_interest(
    i_type = types,
    i = rates,
    m = m
  )

  expect_equal(
    pv_flow(
      cf = cf,
      i = rates,
      i_type = types,
      m = m,
      t = t
    ),
    sum(cf / (1 + i_eff)^t)
  )
})


test_that("pv_flow converts calendar dates using the selected day count", {
  cf <- c(-1000, 1100)
  dates <- as.Date(c("2026-01-01", "2027-01-01"))

  expect_equal(
    pv_flow(
      cf = cf,
      i = 0.10,
      date = dates,
      day_count = "act/365"
    ),
    0,
    tolerance = 1e-10
  )

  t_360 <- c(0, 365 / 360)
  expected_360 <- sum(cf / 1.10^t_360)

  expect_equal(
    pv_flow(
      cf = cf,
      i = 0.10,
      date = dates,
      day_count = "act/360"
    ),
    expected_360
  )
})


test_that("fv_flow agrees with direct accumulation under a constant rate", {
  cf <- c(100, 150, 200)
  t <- c(0, 1, 2)
  i <- 0.08
  T <- max(t)

  expected <- sum(cf * (1 + i)^(T - t))

  expect_equal(
    fv_flow(cf = cf, i = i, t = t),
    expected
  )
})


test_that("constant-rate pv_flow and fv_flow satisfy the focal-date identity", {
  cf <- c(-1000, 300, 400, 500)
  t <- c(0, 1, 2, 3)
  i <- 0.08
  T <- max(t)

  pv <- pv_flow(cf = cf, i = i, t = t)
  fv <- fv_flow(cf = cf, i = i, t = t)

  expect_equal(
    fv,
    pv * (1 + i)^T,
    tolerance = 1e-10
  )
})


test_that("fv_flow uses the horizon spot rate consistently", {
  cf <- c(100, 150, 200)
  t <- c(1, 2, 3)
  spot <- c(0.05, 0.055, 0.06)
  T <- max(t)
  i_T <- spot[t == T][1]

  expected <- sum(
    cf * ((1 + i_T)^T / (1 + spot)^t)
  )

  expect_equal(
    fv_flow(
      cf = cf,
      i = spot,
      t = t
    ),
    expected
  )
})


test_that("fv_flow rejects inconsistent spot rates at a shared horizon", {
  expect_error(
    fv_flow(
      cf = c(100, 200, 300),
      i = c(0.04, 0.05, 0.06),
      t = c(1, 2, 2)
    ),
    "must agree"
  )
})


test_that("general cash-flow valuation validates timing and lengths", {
  expect_error(
    pv_flow(cf = c(100, 200), i = 0.05),
    "either `t` or `date`"
  )

  expect_error(
    fv_flow(cf = c(100, 200), i = 0.05),
    "either `t` or `date`"
  )

  expect_error(
    pv_flow(
      cf = c(100, 200),
      i = 0.05,
      t = c(0, 1),
      date = as.Date(c("2026-01-01", "2027-01-01"))
    ),
    "only one"
  )

  expect_error(
    pv_flow(
      cf = c(100, 200, 300),
      i = c(0.04, 0.05),
      t = c(1, 2, 3)
    ),
    "length 1 or the same length"
  )
})


test_that("empty cash flows have zero present and future value", {
  expect_equal(
    pv_flow(cf = numeric(0), i = 0.05, t = numeric(0)),
    0
  )

  expect_equal(
    fv_flow(cf = numeric(0), i = 0.05, t = numeric(0)),
    0
  )
})
