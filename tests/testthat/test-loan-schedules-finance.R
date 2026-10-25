test_that("amort_schedule level payment matches the actuarial annuity formula", {
  principal <- 100000
  n <- 12
  i <- 0.12
  k <- 12
  i_period <- (1 + i)^(1 / k) - 1
  a_n <- (1 - (1 + i_period)^(-n)) / i_period
  expected_payment <- principal / a_n

  out <- amort_schedule(
    principal = principal,
    n = n,
    i = i,
    k = k,
    timing = "immediate"
  )

  expect_equal(out$payment[1], expected_payment, tolerance = 1e-8)
  expect_equal(out$ob_end[nrow(out)], 0, tolerance = 1e-8)
  expect_equal(sum(out$total_principal), principal, tolerance = 1e-6)
})

test_that("amort_schedule zero-interest loan amortizes linearly", {
  principal <- 12000
  n <- 12

  out <- amort_schedule(
    principal = principal,
    n = n,
    i = 0
  )

  expect_equal(out$payment, rep(principal / n, n), tolerance = 1e-12)
  expect_equal(out$interest, rep(0, n), tolerance = 1e-12)
  expect_equal(out$ob_end[nrow(out)], 0, tolerance = 1e-10)
})

test_that("amort_schedule due payment follows the annuity-due formula", {
  principal <- 50000
  n <- 10
  i <- 0.06
  i_period <- i
  a_n <- (1 - (1 + i_period)^(-n)) / i_period
  expected_payment <- principal / ((1 + i_period) * a_n)

  out <- amort_schedule(
    principal = principal,
    n = n,
    i = i,
    timing = "due"
  )

  expect_equal(out$payment[1], expected_payment, tolerance = 1e-8)
  expect_equal(out$ob_end[nrow(out)], 0, tolerance = 1e-7)
})

test_that("amort_schedule named extra principal is applied in the intended period", {
  out <- amort_schedule(
    principal = 30000,
    n = 12,
    i = 0.08,
    extra_principal = c("4" = 1500),
    adjust = "none"
  )

  expect_equal(out$extra_principal[out$period == 4], 1500)
  expect_equal(sum(out$extra_principal[out$period != 4]), 0)
  expect_equal(out$ob_end[nrow(out)], 0, tolerance = 1e-7)
})

test_that("amort_schedule term adjustment shortens a loan after material prepayment", {
  base <- amort_schedule(
    principal = 100000,
    n = 24,
    i = 0.10,
    k = 12,
    adjust = "none"
  )

  accelerated <- amort_schedule(
    principal = 100000,
    n = 24,
    i = 0.10,
    k = 12,
    extra_principal = c("6" = 25000),
    adjust = "term"
  )

  expect_lt(nrow(accelerated), nrow(base))
  expect_equal(accelerated$ob_end[nrow(accelerated)], 0, tolerance = 1e-7)
  expect_equal(sum(accelerated$total_principal), 100000, tolerance = 1e-5)
})

test_that("amort_schedule payment adjustment preserves the contractual horizon when not prepaid", {
  n <- 24
  out <- amort_schedule(
    principal = 100000,
    n = n,
    i = 0.10,
    k = 12,
    extra_principal = c("6" = 5000),
    adjust = "payment"
  )

  expect_equal(nrow(out), n)
  expect_equal(out$ob_end[nrow(out)], 0, tolerance = 1e-7)
  expect_lt(out$payment[7], out$payment[1])
})

test_that("amort_schedule rejects negative extra principal", {
  expect_error(
    amort_schedule(
      principal = 10000,
      n = 10,
      i = 0.05,
      extra_principal = -100
    ),
    "nonnegative"
  )
})

test_that("sinking_fund_schedule default deposit accumulates exactly to principal", {
  principal <- 100000
  n <- 10
  j <- 0.04
  expected_deposit <- principal / sum((1 + j)^(0:(n - 1)))

  out <- sinking_fund_schedule(
    principal = principal,
    n = n,
    i = 0.07,
    j = j
  )

  expect_equal(out$sinking_deposit[1], expected_deposit, tolerance = 1e-8)
  expect_equal(out$redemption_from_fund[n], principal, tolerance = 1e-6)
  expect_equal(out$loan_balance_end[n], 0, tolerance = 1e-8)
  expect_equal(out$fund_balance_end[n], 0, tolerance = 1e-8)
})

test_that("sinking_fund_schedule loan interest is constant while principal remains outstanding", {
  principal <- 75000
  i <- 0.06
  n <- 8

  out <- sinking_fund_schedule(
    principal = principal,
    n = n,
    i = i,
    j = 0.04
  )

  expect_equal(out$interest_loan, rep(principal * i, n), tolerance = 1e-12)
  expect_equal(out$loan_balance_start, rep(principal, n), tolerance = 1e-12)
})

test_that("sinking_fund_schedule uses equal deposits when the fund rate is zero", {
  principal <- 60000
  n <- 12

  out <- sinking_fund_schedule(
    principal = principal,
    n = n,
    i = 0.05,
    j = 0
  )

  expect_equal(out$sinking_deposit, rep(principal / n, n), tolerance = 1e-12)
  expect_equal(out$fund_balance_end[n], 0, tolerance = 1e-8)
  expect_equal(out$loan_balance_end[n], 0, tolerance = 1e-8)
})

test_that("sinking_fund_schedule converts annual rates consistently to schedule-period rates", {
  out <- sinking_fund_schedule(
    principal = 50000,
    n = 12,
    i = 0.12,
    j = 0.06,
    k = 12
  )

  expect_equal(
    out$i_loan_period[1],
    (1 + 0.12)^(1 / 12) - 1,
    tolerance = 1e-12
  )
  expect_equal(
    out$i_fund_period[1],
    (1 + 0.06)^(1 / 12) - 1,
    tolerance = 1e-12
  )
})

test_that("sinking_fund_schedule rejects nonpositive deposits", {
  expect_error(
    sinking_fund_schedule(
      principal = 10000,
      n = 10,
      i = 0.05,
      j = 0.03,
      deposit = 0
    ),
    "positive number"
  )
})
