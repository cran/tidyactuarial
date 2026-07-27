test_that("constant-rate level schedule agrees with amort_schedule", {
  general <- amort_schedule_general(
    principal = 100000,
    n = 24,
    i = 0.12,
    i_type = "nominal_interest",
    m = 12,
    k = 12
  )

  existing <- amort_schedule(
    principal = 100000,
    n = 24,
    i = 0.12,
    i_type = "nominal_interest",
    m = 12,
    k = 12
  )

  expect_equal(
    general$payment,
    existing$payment,
    tolerance = 1e-8
  )

  expect_equal(
    general$interest,
    existing$interest,
    tolerance = 1e-8
  )

  expect_equal(
    general$ob_end,
    existing$ob_end,
    tolerance = 1e-8
  )
})

test_that("annuity-due level schedule agrees with amort_schedule", {
  general <- amort_schedule_general(
    principal = 50000,
    n = 12,
    i = 0.06,
    timing = "due"
  )

  existing <- amort_schedule(
    principal = 50000,
    n = 12,
    i = 0.06,
    timing = "due"
  )

  expect_equal(
    general$payment,
    existing$payment,
    tolerance = 1e-8
  )

  expect_equal(
    general$ob_end,
    existing$ob_end,
    tolerance = 1e-8
  )
})

test_that("variable rates produce the expected recursive balances", {
  result <- amort_schedule_general(
    principal = 1000,
    n = 3,
    i = c(0.05, 0.06, 0.07),
    payment = c(400, 400, 400)
  )

  expected_end_1 <- 1000 * 1.05 - 400
  expected_end_2 <- expected_end_1 * 1.06 - 400
  expected_end_3 <- max(expected_end_2 * 1.07 - 400, 0)

  expect_equal(
    result$ob_end,
    c(
      expected_end_1,
      expected_end_2,
      expected_end_3
    )
  )
})

test_that("payment vectors are used period by period", {
  payments <- c(100, 200, 300, 400)

  result <- amort_schedule_general(
    principal = 2000,
    n = 4,
    i = 0,
    payment = payments
  )

  expect_equal(
    result$payment,
    payments
  )

  expect_equal(
    result$ob_end,
    2000 - cumsum(payments)
  )
})

test_that("payment functions receive the current loan state", {
  payment_rule <- function(
      period,
      ob_start,
      i_effective_period
  ) {
    100 + 25 * period
  }

  result <- amort_schedule_general(
    principal = 1000,
    n = 4,
    i = 0,
    payment = payment_rule
  )

  expect_equal(
    result$payment,
    c(125, 150, 175, 200)
  )

  expect_true(
    all(result$payment_source == "function")
  )
})

test_that("extra principal is applied without overpaying", {
  result <- amort_schedule_general(
    principal = 1000,
    n = 5,
    i = 0,
    payment = 200,
    extra_principal = c(0, 500, 0, 0, 0)
  )

  expect_equal(
    result$ob_end,
    c(800, 100, 0)
  )

  expect_equal(
    result$extra_principal,
    c(0, 500, 0)
  )
})

test_that("automatically calculated payment amortizes variable-rate loans", {
  result <- amort_schedule_general(
    principal = 10000,
    n = 5,
    i = c(0.04, 0.05, 0.06, 0.07, 0.08),
    payment = NULL
  )

  expect_equal(
    result$ob_end[[nrow(result)]],
    0,
    tolerance = 1e-7
  )

  expect_true(
    all(result$payment_source == "calculated_level")
  )
})

test_that("summary output contains essential actuarial measures", {
  result <- amort_schedule_general(
    principal = 10000,
    n = 12,
    i = 0.06,
    output = "summary"
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)

  expect_named(
    result,
    c(
      "principal",
      "periods_realized",
      "total_interest",
      "total_paid",
      "ending_balance",
      "negative_amortization"
    )
  )

  expect_equal(
    result$ending_balance,
    0,
    tolerance = 1e-7
  )
})

test_that("negative amortization is identified", {
  result <- amort_schedule_general(
    principal = 1000,
    n = 2,
    i = 0.20,
    payment = 50
  )

  expect_true(
    all(result$negative_amortization)
  )

  summary_result <- amort_schedule_general(
    principal = 1000,
    n = 2,
    i = 0.20,
    payment = 50,
    output = "summary"
  )

  expect_true(
    summary_result$negative_amortization
  )
})

test_that("zero-interest schedules are handled exactly", {
  result <- amort_schedule_general(
    principal = 1200,
    n = 12,
    i = 0
  )

  expect_equal(
    result$payment,
    rep(100, 12)
  )

  expect_equal(
    result$interest,
    rep(0, 12)
  )

  expect_equal(
    result$ob_end[[12]],
    0
  )
})

test_that("invalid period-pattern lengths are rejected", {
  expect_error(
    amort_schedule_general(
      principal = 1000,
      n = 3,
      i = c(0.04, 0.05),
      payment = 400
    ),
    "length 1 or `n`"
  )

  expect_error(
    amort_schedule_general(
      principal = 1000,
      n = 3,
      i = 0.05,
      payment = c(200, 300)
    ),
    "length 1 or `n`"
  )
})

test_that("invalid payment-function values are rejected", {
  expect_error(
    amort_schedule_general(
      principal = 1000,
      n = 3,
      i = 0.05,
      payment = function(
        period,
        ob_start,
        i_effective_period
      ) {
        -100
      }
    ),
    "nonnegative"
  )
})
