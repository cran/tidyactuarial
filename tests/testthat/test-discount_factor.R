test_that("discount factor is the reciprocal of accumulation factor", {
  expect_equal(
    discount_factor(
      t = 3,
      i = 0.05
    ),
    1 / 1.05^3
  )
})

test_that("a scalar rate is recycled over several times", {
  times <- c(1, 2, 5)

  expect_equal(
    discount_factor(
      t = times,
      i = 0.07
    ),
    1.07^(-times)
  )
})

test_that("constant-rate results agree with discount_factor_spot", {
  times <- c(1, 2, 4)
  rates <- c(0.04, 0.05, 0.06)

  expect_equal(
    discount_factor(
      t = times,
      i = rates
    ),
    discount_factor_spot(
      t = times,
      i = rates
    )
  )
})

test_that("an accumulation function is evaluated by inverse interval ratios", {
  a_fun <- function(t) {
    exp(0.03 * t + 0.002 * t^2)
  }

  expect_equal(
    discount_factor(
      s = 2,
      t = 5,
      a = a_fun
    ),
    a_fun(2) / a_fun(5)
  )
})

test_that("a time-varying force agrees with its analytic integral", {
  delta_fun <- function(t) {
    0.03 + 0.004 * t
  }

  s <- c(0, 2)
  t <- c(3, 5)

  expected_integral <-
    0.03 * (t - s) +
    0.002 * (t^2 - s^2)

  expect_equal(
    discount_factor(
      s = s,
      t = t,
      delta = delta_fun
    ),
    exp(-expected_integral),
    tolerance = 1e-8
  )
})

test_that("zero-length intervals have discount factor one", {
  expect_equal(
    discount_factor(
      s = c(0, 2, 5),
      t = c(0, 2, 5),
      delta = function(t) 0.04
    ),
    c(1, 1, 1)
  )
})

test_that("tidy output includes accumulation and discount factors", {
  result <- discount_factor(
    s = c(0, 1),
    t = c(2, 4),
    i = 0.06,
    tidy = TRUE
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2L)

  expect_true(
    all(
      c(
        "s",
        "t",
        "elapsed_time",
        "model",
        "accumulation_factor",
        "discount_factor"
      ) %in% names(result)
    )
  )

  expect_equal(
    result$discount_factor,
    1 / result$accumulation_factor
  )
})

test_that("validation is inherited from accumulation_factor", {
  expect_error(
    discount_factor(t = 2),
    "exactly one"
  )

  expect_error(
    discount_factor(
      s = 3,
      t = 2,
      i = 0.05
    ),
    "greater than or equal"
  )

  expect_error(
    discount_factor(
      t = 2,
      i = 0.05,
      delta = function(t) 0.05
    ),
    "exactly one"
  )
})

test_that("missing times propagate", {
  expect_equal(
    discount_factor(
      t = c(1, NA, 3),
      i = 0.05
    ),
    c(1 / 1.05, NA_real_, 1 / 1.05^3)
  )
})
