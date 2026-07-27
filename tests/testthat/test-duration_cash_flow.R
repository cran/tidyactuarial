test_that("duration_cash_flow calculates Macaulay duration correctly", {

  cf <- c(100, 100, 1100)
  t <- c(1, 2, 3)
  rate <- 0.05

  present_values <- cf * (1 + rate)^(-t)

  expected <- sum(t * present_values) /
    sum(present_values)

  result <- duration_cash_flow(
    cf = cf,
    t = t,
    rate = rate,
    duration = "macaulay"
  )

  expect_equal(result, expected)
})


test_that("modified duration is Macaulay duration divided by one plus rate", {

  macaulay <- duration_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = 0.05,
    duration = "macaulay"
  )

  modified <- duration_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = 0.05,
    duration = "modified"
  )

  expect_equal(
    modified,
    macaulay / 1.05
  )
})


test_that("valuation_time excludes earlier and simultaneous cash flows", {

  result <- duration_cash_flow(
    cf = c(-1000, 100, 100, 1100),
    t = c(0, 1, 2, 3),
    rate = 0.05,
    valuation_time = 1,
    duration = "macaulay"
  )

  remaining_times <- c(1, 2)
  present_values <- c(100, 1100) *
    1.05^(-remaining_times)

  expected <- sum(
    remaining_times * present_values
  ) / sum(present_values)

  expect_equal(result, expected)
})


test_that("audit output has the expected structure", {

  result <- duration_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
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
      "duration_contribution"
    )
  )

  expect_equal(nrow(result), 3L)

  expect_equal(
    sum(result$duration_contribution),
    duration_cash_flow(
      cf = c(100, 100, 1100),
      t = c(1, 2, 3),
      rate = 0.05
    )
  )
})


test_that("modified audit contributions sum to modified duration", {

  audit <- duration_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = 0.05,
    duration = "modified",
    output = "audit"
  )

  value <- duration_cash_flow(
    cf = c(100, 100, 1100),
    t = c(1, 2, 3),
    rate = 0.05,
    duration = "modified"
  )

  expect_equal(
    sum(audit$duration_contribution),
    value
  )
})


test_that("cash flows are ordered by time in audit output", {

  result <- duration_cash_flow(
    cf = c(1100, 100, 100),
    t = c(3, 1, 2),
    rate = 0.05,
    output = "audit"
  )

  expect_equal(
    result$time,
    c(1, 2, 3)
  )
})


test_that("a single future payment has duration equal to its remaining time", {

  macaulay <- duration_cash_flow(
    cf = 1000,
    t = 7,
    rate = 0.08,
    valuation_time = 2,
    duration = "macaulay"
  )

  modified <- duration_cash_flow(
    cf = 1000,
    t = 7,
    rate = 0.08,
    valuation_time = 2,
    duration = "modified"
  )

  expect_equal(macaulay, 5)
  expect_equal(modified, 5 / 1.08)
})


test_that("mixed-sign cash flows are accepted when present value is nonzero", {

  result <- duration_cash_flow(
    cf = c(-500, 800),
    t = c(1, 3),
    rate = 0.05
  )

  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(is.finite(result))
})


test_that("duration_cash_flow validates cf and t", {

  expect_error(
    duration_cash_flow(
      cf = c("100", "200"),
      t = c(1, 2),
      rate = 0.05
    ),
    "`cf` and `t` must be numeric vectors",
    fixed = TRUE
  )

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = 1,
      rate = 0.05
    ),
    "`cf` and `t` must have the same length",
    fixed = TRUE
  )

  expect_error(
    duration_cash_flow(
      cf = numeric(),
      t = numeric(),
      rate = 0.05
    ),
    "`cf` and `t` cannot be empty",
    fixed = TRUE
  )

  expect_error(
    duration_cash_flow(
      cf = c(100, Inf),
      t = c(1, 2),
      rate = 0.05
    ),
    "must contain only finite values",
    fixed = TRUE
  )
})


test_that("duration_cash_flow validates rate and valuation_time", {

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = c(0.04, 0.05)
    ),
    "`rate` must be a finite numeric scalar",
    fixed = TRUE
  )

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = -1
    ),
    "`rate` must be greater than -1",
    fixed = TRUE
  )

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      valuation_time = c(0, 1)
    ),
    "`valuation_time` must be a finite numeric scalar",
    fixed = TRUE
  )
})


test_that("an informative error is returned when no future cash flows remain", {

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      valuation_time = 2
    ),
    "No cash flows occur after `valuation_time` = 2",
    fixed = TRUE
  )
})


test_that("duration is undefined when remaining present value is zero", {

  expect_error(
    duration_cash_flow(
      cf = c(100, -105),
      t = c(1, 2),
      rate = 0.05
    ),
    "present value of the remaining cash flows is zero",
    fixed = TRUE
  )
})


test_that("invalid duration and output options are rejected", {

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      duration = "effective"
    )
  )

  expect_error(
    duration_cash_flow(
      cf = c(100, 200),
      t = c(1, 2),
      rate = 0.05,
      output = "table"
    )
  )
})
