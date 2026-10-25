test_that("immunization is invariant to monetary units", {
  for (scale in c(1, 1e-20, 1e20)) {
    out <- immunize_duration(1000 * scale, 5, c(100, 200) * scale,
                              c(3, 7), 0.05)
    expect_equal(out$w, c(1000 / 1.05^5 / 200, 1000 / 1.05^5 / 400),
                 tolerance = 1e-10)
    expect_equal(out$PV_A / out$PV_L, c(1, 1), tolerance = 1e-10)
    three <- immunize_duration_convexity(1000 * scale, 5,
      rep(100, 3) * scale, c(2, 5, 8), c(6, 30, 72) / 1.05^2, 0.05)
    expect_equal(three$w, c(0, 1000 / 1.05^5 / 100, 0), tolerance = 1e-9)
  }
})

test_that("minimum-norm immunization resolves distinct but similar sensitivities", {
  d <- c(9.99, 10, 10.01, 10.02)
  c <- d * (d + 1) / 1.05^2
  out <- immunize_duration_convexity(1000, 10, rep(1, 4), d, c, 0.05)
  pv <- 1000 / 1.05^10
  # Moment equations for offsets (-1,0,1,2) imply the minimum-norm
  # fractions (0.15,0.55,0.45,-0.15); null vector (-1,3,-3,1).
  expect_equal(out$w / pv, c(0.15, 0.55, 0.45, -0.15), tolerance = 1e-7)
  expect_equal(sum(out$w * c) / pv, 110 / 1.05^2, tolerance = 1e-10)
  expect_equal(sum((out$w / pv) * c(-1, 3, -3, 1)), 0, tolerance = 1e-9)
})

test_that("small asset prices do not make distinct durations degenerate", {
  regular <- immunize_duration(1000, 5, c(100, 200), c(3, 7), 0.05)
  tiny <- immunize_duration(1000, 5, c(100, 200) * 1e-20, c(3, 7), 0.05)
  expect_equal(tiny$w * 1e-20, regular$w, tolerance = 1e-12)
})

test_that("duration matching minimizes units rather than market-value fractions", {
  out <- immunize_duration(1000, 5, c(1, 2, 3), c(2, 5, 8), 0.05)
  pv <- 1000 / 1.05^5
  # A null vector of the two constraints is (3,-3,1). Minimum norm is
  # characterized by orthogonality to it, in addition to both constraints.
  expect_equal(sum(out$w * c(1, 2, 3)), pv, tolerance = 1e-10)
  expect_equal(sum(out$w * c(1, 2, 3) * c(2, 5, 8)), 5 * pv,
               tolerance = 1e-10)
  expect_equal(sum(out$w * c(3, -3, 1)) / pv, 0, tolerance = 1e-12)
})

test_that("genuinely singular and undefined immunization systems are rejected", {
  expect_error(immunize_duration(1000, 5, rep(100, 3), rep(3, 3), 0.05),
               "singular")
  expect_error(immunize_duration_convexity(1000, 5, rep(100, 4),
    rep(3, 4), rep(12, 4), 0.05), "singular")
  for (L in list(c(0, 0), c(1000, -1000))) {
    expect_error(immunize_duration(L, c(5, 5), c(100, 200), c(3, 7), 0.05),
                 "Present value")
    expect_error(immunize_duration_convexity(L, c(5, 5), rep(100, 3),
      c(2, 5, 8), c(6, 30, 72), 0.05), "Present value")
  }
  expect_error(immunize_duration(1e308, 1000, c(100, 200), c(3, 7), -0.5),
               "non-finite")
})

test_that("revalued cash flows satisfy duration and convexity matching", {
  i <- 0.05
  maturities <- c(2, 5, 9, 12)
  prices <- 100 / (1 + i)^maturities
  L <- c(300, 700); times <- c(4, 8)
  out <- immunize_duration_convexity(L, times, prices, maturities,
    maturities * (maturities + 1) / (1 + i)^2, i)
  gap <- function(rate) sum(out$w * 100 / (1 + rate)^maturities) -
    sum(L / (1 + rate)^times)
  h <- 1e-4
  expect_equal(gap(i), 0, tolerance = 1e-9)
  expect_lt(abs((gap(i + h) - gap(i - h)) / (2 * h)), 0.001)
  expect_lt(abs((gap(i + h) - 2 * gap(i) + gap(i - h)) / h^2), 0.1)

  duration <- immunize_duration(1000, 5, 100 / (1 + i)^c(3, 7), c(3, 7), i)
  duration_gap <- function(rate) sum(duration$w * 100 / (1 + rate)^c(3, 7)) -
    1000 / (1 + rate)^5
  expect_lt(abs((duration_gap(i + h) - duration_gap(i - h)) / (2 * h)), 0.001)
  expect_gt(duration_gap(i + h), 0)
  expect_gt(duration_gap(i - h), 0)
})

test_that("annual-effective convexity is obtained from coupon sensitivity by the chain rule", {
  i <- 0.05; k <- 2
  b <- bond_convexity(face = 100, c = 0.04, n = 7, k = k, y = i)
  d <- bond_duration(face = 100, c = 0.04, n = 7, k = k, y = i)
  cf <- bond_cash_flows(face = 100, c = 0.04, n = 7, k = k)
  j <- (1 + i)^(1 / k) - 1
  dj <- (1 + j) / (k * (1 + i))
  d2j <- (1 + j) * (1 - k) / (k^2 * (1 + i)^2)
  annual <- b$discrete_convexity_periods * dj^2 -
    d$modified_duration_periods_j * d2j
  expect_equal(annual, convexity_cash_flow(cf$cf, cf$t, i), tolerance = 1e-10)
  expect_gt(abs(annual - b$discrete_convexity_years), 0.1)
})
