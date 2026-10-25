test_that("future_value agrees with the actuarial accumulation formula", {
  expect_equal(
    future_value(C = 1000, i = 0.08, t = 3),
    1000 * 1.08^3
  )
})


test_that("present_value agrees with the actuarial discount formula", {
  expect_equal(
    present_value(C = 1000, i = 0.08, t = 3),
    1000 / 1.08^3
  )
})


test_that("future_value and present_value are inverse operations", {
  capitals <- c(1000, 2500, 5000)
  rates <- c(0.03, 0.06, 0.10)
  times <- c(1, 4, 7)

  accumulated <- future_value(
    C = capitals,
    i = rates,
    t = times
  )

  recovered <- present_value(
    C = accumulated,
    i = rates,
    t = times
  )

  expect_equal(recovered, capitals, tolerance = 1e-10)
})


test_that("single-payment valuation respects nominal interest conversion", {
  j12 <- 0.12
  i_eff <- (1 + j12 / 12)^12 - 1

  expect_equal(
    future_value(
      C = 2000,
      i = j12,
      i_type = "nominal_interest",
      m = 12,
      t = 5
    ),
    2000 * (1 + i_eff)^5
  )

  expect_equal(
    present_value(
      C = 2000,
      i = j12,
      i_type = "nominal_interest",
      m = 12,
      t = 5
    ),
    2000 / (1 + i_eff)^5
  )
})


test_that("single-payment valuation respects a force of interest", {
  delta <- 0.05

  expect_equal(
    future_value(
      C = 1500,
      i = delta,
      i_type = "force",
      t = 4
    ),
    1500 * exp(delta * 4)
  )

  expect_equal(
    present_value(
      C = 1500,
      i = delta,
      i_type = "force",
      t = 4
    ),
    1500 * exp(-delta * 4)
  )
})


test_that("single-payment functions vectorize across rate conventions", {
  C <- c(1000, 1000, 1000, 1000)
  i <- c(0.06, 0.06, 0.06, log(1.06))
  i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )
  m <- c(1, 12, 4, 1)
  t <- c(1, 2, 3, 4)

  i_eff <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  expect_equal(
    future_value(C = C, i = i, i_type = i_type, m = m, t = t),
    C * (1 + i_eff)^t
  )

  expect_equal(
    present_value(C = C, i = i, i_type = i_type, m = m, t = t),
    C / (1 + i_eff)^t
  )
})


test_that("tidy outputs expose actuarially relevant intermediate quantities", {
  fv <- future_value(
    C = c(1000, 2000),
    i = 0.05,
    t = c(2, 3),
    tidy = TRUE
  )

  pv <- present_value(
    C = c(1000, 2000),
    i = 0.05,
    t = c(2, 3),
    tidy = TRUE
  )

  expect_s3_class(fv, "tbl_df")
  expect_s3_class(pv, "tbl_df")

  expect_named(
    fv,
    c(
      "C", "t", "i_input", "i_type", "m", "i_effective",
      "accumulation_factor", "future_value"
    )
  )

  expect_named(
    pv,
    c(
      "C", "t", "i_input", "i_type", "m", "i_effective",
      "v", "present_value"
    )
  )

  expect_equal(fv$future_value, c(1000, 2000) * 1.05^c(2, 3))
  expect_equal(pv$present_value, c(1000, 2000) / 1.05^c(2, 3))
})


test_that("missing values are propagated in single-payment valuation", {
  expect_equal(
    future_value(
      C = c(1000, NA, 3000),
      i = 0.05,
      t = c(1, 2, NA)
    ),
    c(1050, NA_real_, NA_real_)
  )

  expect_equal(
    present_value(
      C = c(1000, NA, 3000),
      i = 0.05,
      t = c(1, 2, NA)
    ),
    c(1000 / 1.05, NA_real_, NA_real_)
  )
})


test_that("single-payment valuation rejects invalid actuarial inputs", {
  expect_error(
    future_value(C = 1000, i = 0.05, t = -1),
    "greater than or equal to 0"
  )

  expect_error(
    present_value(C = 1000, i = 0.05, t = -1),
    "greater than or equal to 0"
  )

  expect_error(
    future_value(
      C = c(1000, 2000),
      i = c(0.04, 0.05, 0.06),
      t = c(1, 2)
    ),
    "length 1 or a common length"
  )

  expect_error(
    present_value(
      C = 1000,
      i = 0.05,
      t = 2,
      tidy = NA
    ),
    "logical scalar"
  )
})
