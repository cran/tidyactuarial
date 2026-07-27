test_that("constant effective interest gives the expected accumulation factor", {
  expect_equal(
    accumulation_factor(
      t = 3,
      i = 0.05
    ),
    1.05^3
  )
})

test_that("a scalar rate is recycled over several times", {
  expect_equal(
    accumulation_factor(
      t = c(1, 2, 5),
      i = 0.07
    ),
    1.07^c(1, 2, 5)
  )
})

test_that("conventional rate specifications reuse standardize_interest", {
  expected_rate <- (1 + 0.12 / 12)^12 - 1

  expect_equal(
    accumulation_factor(
      t = 5,
      i = 0.12,
      i_type = "nominal_interest",
      m = 12
    ),
    (1 + expected_rate)^5
  )
})

test_that("an accumulation function is evaluated by interval ratios", {
  a_fun <- function(t) {
    exp(0.03 * t + 0.002 * t^2)
  }

  expect_equal(
    accumulation_factor(
      s = 2,
      t = 5,
      a = a_fun
    ),
    a_fun(5) / a_fun(2)
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
    accumulation_factor(
      s = s,
      t = t,
      delta = delta_fun
    ),
    exp(expected_integral),
    tolerance = 1e-8
  )
})

test_that("a constant force agrees with the equivalent effective rate", {
  delta_constant <- function(t) {
    log(1.05)
  }

  expect_equal(
    accumulation_factor(
      t = c(1, 3, 7),
      delta = delta_constant
    ),
    accumulation_factor(
      t = c(1, 3, 7),
      i = 0.05
    ),
    tolerance = 1e-8
  )
})

test_that("zero-length intervals have accumulation factor one", {
  expect_equal(
    accumulation_factor(
      s = c(0, 2, 5),
      t = c(0, 2, 5),
      delta = function(t) 0.04
    ),
    c(1, 1, 1)
  )
})

test_that("tidy output contains the actuarially relevant quantities", {
  result <- accumulation_factor(
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
        "i_effective",
        "accumulation_factor"
      ) %in% names(result)
    )
  )
  expect_equal(
    result$accumulation_factor,
    1.06^c(2, 3)
  )
})

test_that("the function requires exactly one interest model", {
  expect_error(
    accumulation_factor(t = 2),
    "exactly one"
  )

  expect_error(
    accumulation_factor(
      t = 2,
      i = 0.05,
      delta = function(t) 0.05
    ),
    "exactly one"
  )
})

test_that("incompatible vector lengths are rejected", {
  expect_error(
    accumulation_factor(
      s = c(0, 1),
      t = c(1, 2, 3),
      i = 0.05
    ),
    "length 1 or a common length"
  )
})

test_that("invalid time order is rejected", {
  expect_error(
    accumulation_factor(
      s = 3,
      t = 2,
      i = 0.05
    ),
    "greater than or equal"
  )
})

test_that("invalid accumulation-function values are rejected", {
  expect_error(
    accumulation_factor(
      t = 2,
      a = function(t) -1
    ),
    "positive finite"
  )
})

test_that("missing times propagate", {
  expect_equal(
    accumulation_factor(
      t = c(1, NA, 3),
      i = 0.05
    ),
    c(1.05, NA_real_, 1.05^3)
  )
})
