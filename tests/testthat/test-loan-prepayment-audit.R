test_that("extras and recasting satisfy dated value and balance identities", {
  for (rate in c(-.02, 0, .08)) {
    for (timing in c("immediate", "due")) {
      for (adjust in c("none", "term", "payment")) {
        out <- amort_schedule(10000, 12, i = rate, timing = timing,
          extra_principal = c("3" = 700, "6" = 300), adjust = adjust)
        b <- as.integer(timing == "due")
        cash_times <- out$period - b
        terminal <- tail(out$period, 1)
        # Independent equation of value, including any outstanding debt.
        value <- sum(out$cashflow * (1 + rate)^(-cash_times)) +
          tail(out$ob_end, 1) * (1 + rate)^(-terminal)
        expect_equal(value, 10000, tolerance = 1e-11)
        expect_equal(out$ob_end, out$ob_start + out$interest - out$cashflow,
          tolerance = 1e-11)
        expect_equal(out$cashflow, out$payment + out$extra_principal)
        expect_true(all(out$payment >= 0 & out$extra_principal >= 0 & out$ob_end >= 0))
        if (adjust == "payment") {
          expect_equal(nrow(out), 12L)
          # Each next regular payment values the remaining balance, before
          # incorporating any future, unanticipated extra payment.
          for (j in seq_len(nrow(out) - 1L)) {
            times <- seq_len(12L - j) - b
            expect_equal(out$payment[j + 1L] * sum((1 + rate)^(-times)),
              out$ob_end[j], tolerance = 1e-11)
          }
        }
      }
    }
  }
})

test_that("payoff stops the schedule and preserves each local allocation convention", {
  for (timing in c("immediate", "due")) {
    fixed <- amort_schedule(100, 5, i = .1, payment = 20,
      extra_principal = 1000, timing = timing)
    general <- amort_schedule_general(100, 5, i = .1, payment = 20,
      extra_principal = 1000, timing = timing)
    due <- if (timing == "immediate") 110 else 100
    expect_equal(nrow(fixed), 1L)
    expect_equal(nrow(general), 1L)
    expect_equal(fixed$cashflow, due)
    expect_equal(general$cashflow, due)
    expect_equal(fixed$ob_end, 0)
    expect_equal(general$ob_end, 0)
    # Fixed schedule preserves requested extras first on the payoff row;
    # the general schedule preserves the requested regular payment first.
    expect_equal(fixed$payment, 0)
    expect_equal(general$payment, 20)
  }
  calls <- integer()
  out <- amort_schedule_general(100, 6, i = 0,
    payment = function(period, ob_start, i_effective_period) {
      calls <<- c(calls, period)
      200
    })
  expect_equal(calls, 1L)
  expect_equal(out$cashflow, 100)
})

test_that("named extras aggregate exact periods without truncating fractions", {
  duplicated <- amort_schedule(10000, 12, i = .04,
    extra_principal = stats::setNames(c(100, 200), c("3", "3")))
  dense <- amort_schedule(10000, 12, i = .04,
    extra_principal = c(0, 0, 300, rep(0, 9)))
  expect_equal(duplicated, dense)
  for (name in c("2.9", "1.5", "2oops", "", "Inf")) {
    expect_error(amort_schedule(10000, 12, i = .04,
      extra_principal = stats::setNames(100, name)), "Names of `extra_principal`")
  }
})

test_that("negative amortization in advance payments is identified", {
  out <- amort_schedule_general(100, 2, i = .1, timing = "due", payment = 1)
  expect_equal(out$ob_end, c(108.9, 118.69), tolerance = 1e-12)
  expect_true(all(out$negative_amortization))
  summary <- amort_schedule_general(100, 2, i = .1, timing = "due", payment = 1,
    output = "summary")
  expect_true(summary$negative_amortization)
  expect_equal(summary$ending_balance,
    100 + summary$total_interest - summary$total_paid, tolerance = 1e-12)
  covered <- amort_schedule_general(100, 2, i = .1, timing = "due", payment = 1,
    extra_principal = 10)
  expect_false(any(covered$negative_amortization))
})

test_that("variable rates, summary and callback payments retain the same ledger", {
  rates <- c(.03, -.01, .08, .04)
  args <- list(principal = 1000, n = 4, i = rates,
    payment = c(100, 200, 150, 800), extra_principal = c(0, 50, 0, 0))
  for (timing in c("immediate", "due")) {
    out <- do.call(amort_schedule_general, c(args, list(timing = timing)))
    summary <- do.call(amort_schedule_general, c(args,
      list(timing = timing, output = "summary")))
    discount <- 1 / cumprod(1 + rates)
    payment_discount <- if (timing == "immediate") discount else c(1, head(discount, -1))
    expect_equal(sum(out$cashflow * payment_discount[out$period]) +
      tail(out$ob_end, 1) * discount[tail(out$period, 1)], 1000, tolerance = 1e-12)
    expect_equal(summary$total_paid, sum(out$cashflow))
    expect_equal(summary$total_interest, sum(out$interest))
    expect_equal(summary$ending_balance, tail(out$ob_end, 1))
    expect_equal(out$ob_end, out$ob_start + out$interest - out$cashflow,
      tolerance = 1e-12)
  }
})

test_that("balance tolerances must be finite", {
  expect_error(amort_schedule(1000, 3, i = .04, tol = Inf), "finite")
  expect_error(sinking_fund_schedule(1000, 3, i = .04, j = .02, tol = Inf), "finite")
})

test_that("representable calculated payments survive overflowing discount sums", {
  # Decimal 90: 2^1000*(1-exp(-2))/(exp(800)-1), payments in arrears.
  expected <- 3.39826998028498449755325673227064986474e-47
  for (fun in list(amort_schedule, amort_schedule_general)) {
    out <- fun(2^1000, 400, i = -2, i_type = "force", tol = 1e-100)
    expect_equal(out$payment[1] / expected, 1, tolerance = 1e-12)
    expect_true(all(is.finite(out$ob_end)))
  }
})

test_that("nonrepresentable payments and monetary ledgers fail explicitly", {
  for (fun in list(amort_schedule, amort_schedule_general)) {
    expect_error(fun(1e-300, 400, i = -2, i_type = "force", tol = 1e-310),
      "Calculated level payment")
    expect_error(fun(1e308, 3, i = 10, payment = 1), "monetary range")
  }
  expect_error(sinking_fund_schedule(100, 3, i = .04, j = 10, deposit = 1e308),
    "monetary range")
})

test_that("a supplied initial payment is recast at the next contractual step", {
  out <- amort_schedule(1000, 4, i = .1, payment = 100, adjust = "payment")
  expect_equal(out$payment[1], 100)
  expect_equal(out$ob_end[1], 1000)
  expect_equal(out$payment[2] * sum(1.1^(-(1:3))), 1000, tolerance = 1e-12)
  expect_equal(tail(out$ob_end, 1), 0, tolerance = 1e-8)
})

test_that("custom sinking deposits retain underfunding and surplus transparently", {
  for (deposit in c(25, 40, 75)) {
    if (deposit == 40) {
      out <- expect_no_warning(sinking_fund_schedule(100, 2, i = .1, j = .5,
        deposit = deposit))
    } else {
      expect_warning(out <- sinking_fund_schedule(100, 2, i = .1, j = .5,
        deposit = deposit), "Final mismatch")
    }
    accumulated <- 2.5 * deposit
    expect_equal(tail(out$redemption_from_fund, 1), min(accumulated, 100))
    expect_equal(tail(out$loan_balance_end, 1), max(100 - accumulated, 0))
    expect_equal(tail(out$fund_balance_end, 1), max(accumulated - 100, 0),
      tolerance = 1e-12)
    expect_equal(out$total_cashflow_borrower, rep(10 + deposit, 2))
    expect_equal(out$fund_balance_end_before_redemption,
      out$fund_balance_start + out$interest_fund + out$sinking_deposit)
  }
})
