test_that("cash-flow sensitivities agree with price derivatives and contributions", {
  times <- c(0, 0.5, 1.5, 4, 7)
  for (payments in list(c(999, 40, 60, 80, 500), c(999, -20, 60, 80, 500))) {
    for (rate in c(-0.03, 0, 0.07)) {
      for (valuation in c(0, 1)) {
        active <- times > valuation
        dt <- times[active] - valuation
        cash <- payments[active]
        p <- sum(cash / (1 + rate)^dt)
        expected_d <- sum(dt * cash / (1 + rate)^(dt + 1)) / p
        expected_c <- sum(dt * (dt + 1) * cash / (1 + rate)^(dt + 2)) / p
        d <- duration_cash_flow(payments, times, rate, valuation, duration = "modified")
        c <- convexity_cash_flow(payments, times, rate, valuation)
        expect_equal(d, expected_d, tolerance = 2e-13)
        expect_equal(c, expected_c, tolerance = 2e-13)
        price <- function(x) sum(cash / (1 + x)^dt)
        h <- 1e-4
        expect_equal(-(price(rate + h) - price(rate - h)) / (2 * h * p),
                     d, tolerance = 2e-6)
        expect_equal((price(rate + h) - 2 * p + price(rate - h)) / (h^2 * p),
                     c, tolerance = 2e-6)
        ad <- duration_cash_flow(payments, times, rate, valuation,
                                duration = "modified", output = "audit")
        ac <- convexity_cash_flow(payments, times, rate, valuation, output = "audit")
        expect_equal(sum(ad$duration_contribution), d, tolerance = 2e-13)
        expect_equal(sum(ac$convexity_contribution), c, tolerance = 2e-13)
        expect_true(all(ad$time > valuation))
      }
    }
  }
})

test_that("cash-flow moments avoid overflowing monetary numerators", {
  for (scale in c(1e-250, 1, 1e308)) {
    expect_equal(duration_cash_flow(scale, 2, 0), 2, tolerance = 3e-13)
    expect_equal(convexity_cash_flow(scale, 2, 0), 6, tolerance = 3e-13)
  }
  audit <- duration_cash_flow(1e308, 2, 0.5, duration = "modified", output = "audit")
  expect_equal(audit$duration_contribution, 2 / 1.5, tolerance = 3e-13)
  # The discount underflows by itself but its product with the amount is finite.
  expect_equal(duration_cash_flow(1e300, 2000, 1), 2000, tolerance = 3e-13)
  expect_equal(convexity_cash_flow(1e300, 2000, 1), 2000 * 2001 / 4,
               tolerance = 3e-13)
  # A weighted moment can survive even if PV/total_PV underflows first.
  cf <- c(1e200, 1e-200)
  times <- c(1e-300, 1e200)
  expect_equal(duration_cash_flow(cf, times, 0) / 1e-200, 1, tolerance = 3e-13)
  expect_equal(convexity_cash_flow(cf, times, 0), 1, tolerance = 3e-13)
  expect_equal(duration_cash_flow(-1e308, 2, 0), 2, tolerance = 3e-13)
  expect_equal(convexity_cash_flow(-1e308, 2, 0), 6, tolerance = 3e-13)
  expect_error(duration_cash_flow(c(1, -1 + 1e-10), c(1, 2), 0), "close to zero")
  expect_error(convexity_cash_flow(c(1, -1 + 1e-10), c(1, 2), 0), "close to zero")
})

test_that("bond sensitivities retain currency scale and discount range", {
  for (face in c(1e-250, 100, 1e308)) {
    d <- bond_duration(face, 0, 2, y = 0.05)
    c <- bond_convexity(face, 0, 2, y = 0.05)
    expect_equal(d$macaulay_duration_years, 2, tolerance = 3e-13)
    expect_equal(d$modified_duration_years_i, 2 / 1.05, tolerance = 3e-13)
    expect_equal(c$discrete_convexity_years, 6 / 1.05^2, tolerance = 3e-13)
  }
  d <- bond_duration(1e300, 0, 2, y_effective_per_period = 1e200)
  expect_equal(d$price / 1e-100, 1, tolerance = 3e-13)
  expect_equal(d$macaulay_duration_years, 2, tolerance = 3e-13)
  d <- bond_duration(100, 0, 0.01, k = 100, y = -100, y_type = "force")
  expect_equal(d$price, 100 * exp(1), tolerance = 3e-13)
  expect_equal(d$modified_duration_years_i / (0.01 * exp(100)), 1,
               tolerance = 3e-13)
  c <- bond_convexity(100, 0, 0.01, k = 100, y = -100, y_type = "force")
  expect_equal(c$discrete_convexity_periods, 2 * exp(2), tolerance = 3e-13)
  expect_no_warning(expect_error(bond_duration(100, 0.05, 2^31, y = 0.05), "n \\* k"))
  expect_no_warning(expect_error(bond_convexity(100, 0.05, 2^31, y = 0.05), "n \\* k"))
  expect_error(bond_duration(0, 0, 2, y = 0.05), "price is zero")
  expect_error(bond_convexity(0, 0, 2, y = 0.05), "price is zero")
})

test_that("bond derivatives distinguish coupon-period, nominal and annual yields", {
  for (k in c(1L, 2L, 4L)) {
    for (j in c(-0.01, 0.025)) {
      d <- bond_duration(1000, 0.06, 7, k, y_effective_per_period = j)
      c <- bond_convexity(1000, 0.06, 7, k, y_effective_per_period = j)
      periods <- c(seq_len(7 * k), 7 * k)
      cash <- c(rep(60 / k, 7 * k), 1000)
      price <- function(x) sum(cash / (1 + x)^periods)
      p <- price(j)
      h <- 1e-5
      expect_equal(-(price(j + h) - price(j - h)) / (2 * h * p),
                   d$modified_duration_periods_j, tolerance = 2e-6)
      expect_equal((price(j + h) - 2 * p + price(j - h)) / (h^2 * p),
                   c$discrete_convexity_periods, tolerance = 2e-6)
      annual <- (1 + j)^k - 1
      times <- periods / k
      pv <- cash / (1 + annual)^times
      annual_c <- sum(times * (times + 1) * pv) / sum(pv) / (1 + annual)^2
      dj <- (1 + j) / (k * (1 + annual))
      d2j <- (1 + j) * (1 - k) / (k^2 * (1 + annual)^2)
      chain_c <- c$discrete_convexity_periods * dj^2 -
        d$modified_duration_periods_j * d2j
      expect_equal(chain_c, annual_c, tolerance = 3e-13)
      expect_equal(c$discrete_convexity_years, c$discrete_convexity_periods / k^2,
                   tolerance = 3e-13)
      expect_equal(d$modified_duration_years_i,
                   sum(times * pv) / sum(pv) / (1 + annual), tolerance = 3e-13)
    }
  }
})

test_that("portfolio sensitivities avoid monetary overflow and underflow", {
  for (scale in c(1e-300, 1, 1e307)) {
    expect_equal(portfolio_duration(P = rep(scale, 2), D = c(100, 200))$D_P,
                 150, tolerance = 3e-13)
    expect_equal(portfolio_convexity(P = rep(scale, 2), C = c(1000, 2000))$C_P,
                 1500, tolerance = 3e-13)
  }
  expect_equal(portfolio_convexity(P = c(1e-300, 1e-300), C = c(1e-100, 2e-100))$C_P /
                 1e-100, 1.5, tolerance = 3e-13)
  expect_equal(portfolio_duration(P = c(1e-300, 1e300), D = c(1e300, 0))$D_P /
                 1e-300, 1, tolerance = 3e-13)
  expect_equal(portfolio_duration(P = c(0, 100), D = c(1e308, 4))$D_P, 4)
  expect_error(portfolio_duration(P = c(1e308, 1e308), D = c(1, 2)), "total value")
  expect_error(portfolio_convexity(P = c(0, 0), C = c(1, 2)), "total value")
})

test_that("portfolio aggregation represents derivatives under a common yield basis", {
  rate <- 0.04
  schedules <- list(list(t = 1:4, cf = c(5, 5, 5, 105)),
                    list(t = c(0.5, 1.5, 5), cf = c(20, 30, 220)))
  values <- vapply(schedules, function(s) sum(s$cf / (1 + rate)^s$t), numeric(1))
  duration <- vapply(schedules, function(s)
    sum(s$t * s$cf / (1 + rate)^(s$t + 1)) /
      sum(s$cf / (1 + rate)^s$t), numeric(1))
  convexity <- vapply(schedules, function(s)
    sum(s$t * (s$t + 1) * s$cf / (1 + rate)^(s$t + 2)) /
      sum(s$cf / (1 + rate)^s$t), numeric(1))
  d <- portfolio_duration(P = values, D = duration)$D_P
  c <- portfolio_convexity(P = values, C = convexity)$C_P
  price <- function(r) sum(vapply(schedules, function(s)
    sum(s$cf / (1 + r)^s$t), numeric(1)))
  h <- 1e-4
  p <- price(rate)
  expect_equal(-(price(rate + h) - price(rate - h)) / (2 * h * p), d,
               tolerance = 2e-6)
  expect_equal((price(rate + h) - 2 * p + price(rate - h)) / (h^2 * p), c,
               tolerance = 2e-6)
  error <- function(h) (price(rate + h) - p) / p - (-d * h + c * h^2 / 2)
  expect_lt(abs(error(0.0005)), 0.2 * abs(error(0.001)))
})
