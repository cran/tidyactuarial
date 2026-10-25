test_that("discount_factor_spot matches the actuarial spot discount formula", {
  t <- c(0, 1, 2, 5)
  i <- c(0.03, 0.04, 0.05, 0.06)

  expected <- (1 + i)^(-t)
  actual <- discount_factor_spot(t = t, i = i)

  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(actual[1], 1)
})


test_that("discount_factor_spot standardizes nominal interest before discounting", {
  t <- c(1, 2, 4)
  j <- c(0.06, 0.08, 0.10)
  m <- 2

  i_eff <- (1 + j / m)^m - 1
  expected <- (1 + i_eff)^(-t)

  actual <- discount_factor_spot(
    t = t,
    i = j,
    i_type = "nominal_interest",
    m = m
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})


test_that("discount_factor_spot standardizes force of interest before discounting", {
  t <- c(1, 2, 3)
  delta <- c(0.03, 0.04, 0.05)

  i_eff <- exp(delta) - 1
  expected <- exp(-delta * t)

  actual <- discount_factor_spot(
    t = t,
    i = delta,
    i_type = "force"
  )

  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(actual, (1 + i_eff)^(-t), tolerance = 1e-12)
})


test_that("discount_factor_spot is vectorized and recycles scalar inputs", {
  t <- c(1, 2, 3, 4)
  i <- 0.05

  expect_equal(
    discount_factor_spot(t = t, i = i),
    (1 + i)^(-t),
    tolerance = 1e-12
  )
})


test_that("discount_factor_spot tidy output exposes the audit trail", {
  out <- discount_factor_spot(
    t = c(1, 2),
    i = c(0.08, 0.10),
    i_type = "nominal_interest",
    m = 2,
    tidy = TRUE
  )

  expected_i <- (1 + c(0.08, 0.10) / 2)^2 - 1

  expect_s3_class(out, "tbl_df")
  expect_named(
    out,
    c("t", "i_input", "i_type", "m", "i_effective", "discount_factor")
  )
  expect_equal(out$i_effective, expected_i, tolerance = 1e-12)
  expect_equal(
    out$discount_factor,
    (1 + expected_i)^(-out$t),
    tolerance = 1e-12
  )
})


test_that("discount_factor_spot propagates missing values", {
  out <- discount_factor_spot(
    t = c(1, 2, NA_real_),
    i = c(0.04, NA_real_, 0.06)
  )

  expect_equal(out[1], 1 / 1.04, tolerance = 1e-12)
  expect_true(is.na(out[2]))
  expect_true(is.na(out[3]))
})


test_that("discount_factor_spot validates actuarially invalid inputs", {
  expect_error(
    discount_factor_spot(t = -1, i = 0.05),
    "greater than or equal to 0"
  )

  expect_error(
    discount_factor_spot(t = 2, i = -1.01),
    "greater than -1"
  )

  expect_error(
    discount_factor_spot(t = c(1, 2), i = c(0.04, 0.05, 0.06)),
    "length 1 or a common length"
  )
})
