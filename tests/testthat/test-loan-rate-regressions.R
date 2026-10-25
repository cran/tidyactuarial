test_that("loan balance tolerance is not a zero-interest threshold", {
  principal <- 1e9
  n <- 24L
  for (delta in c(5e-9, .005)) {
    monetary_tol <- if (delta == .005) .01 else 1e-8
    for (timing in c("immediate", "due")) {
      times <- seq_len(n) - as.integer(timing == "due")
      expected <- principal / sum(exp(-delta * times))
      out <- amort_schedule(principal, n, i = delta, i_type = "force",
        timing = timing, tol = monetary_tol)
      expect_equal(out$payment[1] / expected, 1, tolerance = 1e-13)
      # Reconstruct balances from the original rate specification.
      balance <- principal
      for (row in seq_len(nrow(out))) {
        balance <- if (timing == "immediate") {
          balance * exp(delta) - out$cashflow[row]
        } else {
          (balance - out$cashflow[row]) * exp(delta)
        }
        expect_equal(out$ob_end[row] / principal, balance / principal, tolerance = 1e-13)
      }
      expect_lt(abs(tail(out$ob_end, 1)) / principal, 1e-13)
    }
  }
})

test_that("finite period forces need no representable annual equivalent", {
  delta <- 800
  k <- 1200
  n <- 3L
  expected <- 1000 / sum(exp(-delta * seq_len(n) / k))
  regular <- amort_schedule(1000, n, i = delta, i_type = "force", k = k)
  general <- amort_schedule_general(1000, n, i = delta, i_type = "force", k = k)
  expect_equal(regular$payment[1], expected, tolerance = 1e-12)
  expect_equal(general$payment[1], expected, tolerance = 1e-12)
  expect_equal(regular$ob_end, general$ob_end, tolerance = 1e-11)
})

test_that("sinking-fund deposits retain small fund interest", {
  principal <- 1e9
  n <- 24L
  for (delta in c(5e-9, .005)) {
    expected <- principal / sum(exp(delta * (0:(n - 1L))))
    out <- expect_no_warning(sinking_fund_schedule(principal, n, i = .04,
      j = delta, j_type = "force", tol = if (delta == .005) .01 else 1e-4))
    expect_equal(out$sinking_deposit[1] / expected, 1, tolerance = 1e-13)
    expect_lt(abs(tail(out$loan_balance_end, 1)) / principal, 1e-13)
    expect_lt(abs(tail(out$fund_balance_end, 1)) / principal, 1e-13)
  }
})
