test_that("convexity agrees with its defining formula", {

  cf <- c(100, 100, 1100)
  t <- c(1, 2, 3)
  rate <- 0.05

  present_value <- cf *
    (1 + rate)^(-t)

  expected <- sum(
    t *
      (t + 1) *
      present_value /
      (1 + rate)^2
  ) / sum(present_value)

  result <- convexity_cash_flow(
    cf = cf,
    t = t,
    rate = rate
  )

  expect_equal(result, expected)
  expect_type(result, "double")
  expect_length(result, 1L)
})


test_that("convexity agrees with numerical second price sensitivity", {

  cf <- c(80, 80, 80, 1080)
  t <- c(0.5, 1.5, 2.5, 3.5)
  valuation_time <- 0.25
  rate <- 0.06
  h <- 1e-3

  price_at <- function(yield) {
    remaining_time <- t - valuation_time

    sum(
      cf *
        (1 + yield)^(-remaining_time)
    )
  }

  second_derivative <- (
    -price_at(rate + 2 * h) +
      16 * price_at(rate + h) -
      30 * price_at(rate) +
      16 * price_at(rate - h) -
      price_at(rate - 2 * h)
  ) / (12 * h^2)

  numerical_convexity <-
    second_derivative /
      price_at(rate)

  result <- convexity_cash_flow(
    cf = cf,
    t = t,
    rate = rate,
    valuation_time = valuation_time
  )

  expect_equal(
    result,
    numerical_convexity,
    tolerance = 1e-7
  )
})


test_that("valuation_time uses remaining times and excludes simultaneous flows", {

  result <- convexity_cash_flow(
    cf = c(-1000, 5000, 100, 1100),
    t = c(0, 1, 2, 3),
    rate = 0.05,
    valuation_time = 1
  )

  remaining_time <- c(1, 2)

  present_value <- c(100, 1100) *
    1.05^(-remaining_time)

  expected <- sum(
    remaining_time *
      (remaining_time + 1) *
      present_value /
      1.05^2
  ) / sum(present_value)

  expect_equal(result, expected)
})


test_that("convexity is invariant to a common translation of all times", {

  original <- convexity_cash_flow(
    cf = c(100, 100, 1100),
    t = c(2, 3, 5),
    rate = 0.04,
    valuation_time = 1
  )

  translated <- convexity_cash_flow(
    cf = c(100, 100, 1100),
    t = c(12, 13, 15),
    rate = 0.04,
    valuation_time = 11
  )

  expect_equal(translated, original)
})


test_that("at zero yield convexity is the weighted value of u times u plus one", {

  cf <- c(100, 200, 300)
  t <- c(1, 2, 4)

  result <- convexity_cash_flow(
    cf = cf,
    t = t,
    rate = 0
  )

  expected <- sum(
    t *
      (t + 1) *
      cf
  ) / sum(cf)

  expect_equal(result, expected)
})


test_that("a single payment has the correct convexity", {

  result <- convexity_cash_flow(
    cf = 1000,
    t = 7.5,
    rate = 0.08,
    valuation_time = 2
  )

  remaining_time <- 5.5

  expected <- remaining_time *
    (remaining_time + 1) /
    1.08^2

  expect_equal(result, expected)
})


test_that("second-order duration-convexity approximation is locally accurate", {

  cf <- c(60, 60, 60, 1060)
  t <- c(1, 2, 3, 4)
  rate <- 0.05
  delta_rate <- 1e-4

  price_at <- function(yield) {
    sum(
      cf *
        (1 + yield)^(-t)
    )
  }

  modified_duration <- duration_cash_flow(
    cf = cf,
    t = t,
    rate = rate,
    duration = "modified"
  )

  convexity <- convexity_cash_flow(
    cf = cf,
    t = t,
    rate = rate
  )

  exact_relative_change <- (
    price_at(rate + delta_rate) -
      price_at(rate)
  ) / price_at(rate)

  approximate_relative_change <-
    -modified_duration * delta_rate +
    0.5 * convexity * delta_rate^2

  approximation_error <- abs(
    approximate_relative_change -
      exact_relative_change
  )

  expect_lt(
    approximation_error,
    5e-10
  )
})

test_that("audit output has the documented structure", {

  result <- convexity_cash_flow(
    cf = c(1100, 100, 100),
    t = c(3, 1, 2),
    rate = 0.05,
    output = "audit"
  )

  expect_s3_class(result, "tbl_df")

  expect_named(
    result,
    c(
      "time",
      "remaining_time",
      "cash_flow",
      "discount_factor",
      "present_value",
      "convexity_contribution"
    )
  )

  expect_equal(result$time, c(1, 2, 3))
  expect_equal(nrow(result), 3L)
})


test_that("audit contributions sum to convexity", {

  audit <- convexity_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = 0.05,
    output = "audit"
  )

  value <- convexity_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = 0.05,
    output = "value"
  )

  expect_equal(
    sum(audit$convexity_contribution),
    value
  )
})


test_that("duplicate payment times are retained in stable order", {

  result <- convexity_cash_flow(
    cf = c(300, 100, 200),
    t = c(2, 1, 2),
    rate = 0.03,
    output = "audit"
  )

  expect_equal(result$time, c(1, 2, 2))
  expect_equal(result$cash_flow, c(100, 300, 200))
})


test_that("mixed-sign cash flows are accepted when present value is nonzero", {

  result <- convexity_cash_flow(
    cf = c(-500, 800),
    t = c(1, 3),
    rate = 0.05
  )

  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(is.finite(result))
})


test_that("mixed-sign cash flows may produce negative convexity", {

  result <- convexity_cash_flow(
    cf = c(1000, -900),
    t = c(1, 5),
    rate = 0.05
  )

  expect_lt(result, 0)
})


test_that("valid negative yields greater than minus one are accepted", {

  result <- convexity_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = -0.02
  )

  expect_true(is.finite(result))
})


test_that("small monetary scales are not mistaken for zero present value", {

  result <- convexity_cash_flow(
    cf = c(1e-20, 2e-20),
    t = c(1, 2),
    rate = 0.05
  )

  expect_true(is.finite(result))
})


test_that("cf and t must be finite numeric vectors of equal positive length", {

  expect_error(
    convexity_cash_flow(
      cf = c("100", "200"),
      t = c(1, 2),
      rate = 0.05
    ),
    "`cf` and `t` must be numeric vectors",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = matrix(c(100, 200), ncol = 1),
      t = c(1, 2),
      rate = 0.05
    ),
    "`cf` and `t` must be numeric vectors",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = 1,
      rate = 0.05
    ),
    "`cf` and `t` must have the same length",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = numeric(),
      t = numeric(),
      rate = 0.05
    ),
    "`cf` and `t` cannot be empty",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, Inf),
      t = c(1, 2),
      rate = 0.05
    ),
    "must contain only finite values",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, NA_real_),
      rate = 0.05
    ),
    "must contain only finite values",
    fixed = TRUE
  )
})


test_that("rate and valuation_time must be valid numeric scalars", {

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = c(0.04, 0.05)
    ),
    "`rate` must be a finite numeric scalar",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = NA_real_
    ),
    "`rate` must be a finite numeric scalar",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = -1
    ),
    "`rate` must be greater than -1",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      valuation_time = c(0, 1)
    ),
    "`valuation_time` must be a finite numeric scalar",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      valuation_time = Inf
    ),
    "`valuation_time` must be a finite numeric scalar",
    fixed = TRUE
  )
})


test_that("an informative error is returned when no future cash flows remain", {

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      valuation_time = 2
    ),
    "No cash flows occur after `valuation_time` = 2",
    fixed = TRUE
  )
})


test_that("convexity is undefined for zero or cancelling present value", {

  expect_error(
    convexity_cash_flow(
      cf = c(0, 0),
      t = c(1, 2),
      rate = 0.05
    ),
    "present value of the remaining cash flows is zero",
    fixed = TRUE
  )

  expect_error(
    convexity_cash_flow(
      cf = c(100, -105),
      t = c(1, 2),
      rate = 0.05
    ),
    "present value of the remaining cash flows is zero",
    fixed = TRUE
  )
})


test_that("non-finite discounting is detected explicitly", {

  expect_error(
    convexity_cash_flow(
      cf = 100,
      t = 1e308,
      rate = -0.5
    ),
    "Discounting produced non-finite values",
    fixed = TRUE
  )
})


test_that("non-finite convexity weights are detected explicitly", {

  expect_error(
    convexity_cash_flow(
      cf = 1,
      t = 1e200,
      rate = 1e-310
    ),
    "convexity weights are non-finite",
    fixed = TRUE
  )
})


test_that("invalid output options are rejected", {

  expect_error(
    convexity_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      output = "table"
    )
  )
})
