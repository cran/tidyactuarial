test_that("immunization preserves finite moments and systems across monetary scales", {
  for (scale in c(1e-250, 1, 1e305)) {
    for (prices in list(c(100, 200), c(100, 200, 300))) {
      durations <- if (length(prices) == 2) c(3, 7) else c(2, 5, 8)
      out <- immunize_duration(1000 * scale, 5, prices * scale, durations, 0)
      expect_equal(out$PV_L / (1000 * scale), rep(1, length(prices)), tolerance = 3e-12)
      expect_equal(out$PV_A / out$PV_L, rep(1, length(prices)), tolerance = 3e-12)
      expect_equal(out$D_L, rep(5, length(prices)), tolerance = 3e-12)
      expect_equal(out$D_A, rep(5, length(prices)), tolerance = 3e-12)
      expect_equal(sum(out$w * prices), 1000, tolerance = 3e-12)
      expect_equal(sum(out$w * prices * durations), 5000, tolerance = 3e-12)
      if (length(prices) == 2) {
        expect_equal(out$w, c(5, 2.5), tolerance = 3e-12)
      } else {
        # Null vector of the original two constraints: units, not value fractions.
        expect_equal(sum(out$w * c(3, -3, 1)), 0, tolerance = 3e-12)
      }
    }
    for (durations in list(c(2, 5, 8), c(2, 4, 6, 8))) {
      convexities <- durations * (durations + 1)
      out <- immunize_duration_convexity(1000 * scale, 5,
                 rep(100 * scale, length(durations)), durations, convexities, 0)
      expect_equal(out$PV_A / out$PV_L, rep(1, length(durations)), tolerance = 3e-12)
      expect_equal(out$D_A, rep(5, length(durations)), tolerance = 3e-12)
      expect_equal(out$C_A, rep(30, length(durations)), tolerance = 3e-12)
      expect_equal(sum(out$w * 100), 1000, tolerance = 3e-12)
      expect_equal(sum(out$w * 100 * durations), 5000, tolerance = 3e-12)
      expect_equal(sum(out$w * 100 * convexities), 30000, tolerance = 3e-12)
      if (length(durations) == 4) {
        expect_equal(sum(out$w * c(-1, 3, -3, 1)), 0, tolerance = 3e-12)
      }
    }
  }
})

test_that("liability valuation retains original growth and derivative coordinates", {
  out <- immunize_duration(100, 1e20, c(1, 1), c(0.5e20, 1.5e20), 1e-20)
  expect_equal(out$PV_L, rep(100 * exp(-1), 2), tolerance = 3e-13)
  out <- immunize_duration_convexity(100, 0.01, c(1, 1, 1),
             c(0, 0.01, 0.02), c(0, 1, 4) * exp(200), -100, i_type = "force")
  expect_equal(out$PV_L, rep(100 * exp(1), 3), tolerance = 3e-12)
  expect_equal(out$D_L, rep(0.01, 3), tolerance = 3e-12)
  expect_equal(out$C_L / exp(200), rep(0.01 * 1.01, 3), tolerance = 3e-12)
  expect_equal(out$C_A / out$C_L, rep(1, 3), tolerance = 3e-12)
  expect_no_warning(x <- immunize_duration(100, 5, c(10, 10), c(3, 7), 0.05, m = 2^31))
  expect_equal(x$PV_L, rep(100 / 1.05^5, 2), tolerance = 3e-13)
  x <- immunize_duration(-1000, 5, c(100, 200), c(3, 7), 0)
  expect_equal(x$w, c(-5, -2.5), tolerance = 3e-12)
  expect_equal(x$PV_A, rep(-1000, 2), tolerance = 3e-12)
  x <- immunize_duration(100, 0, c(10, 20), c(0, 2), 0.05)
  expect_equal(x$w, c(10, 0), tolerance = 3e-12)
  expect_equal(x$D_L, c(0, 0))
  x <- immunize_duration_convexity(100, 0, c(10, 20, 30),
             c(0, 2, 5), c(0, 6, 30) / 1.05^2, 0.05)
  expect_equal(x$w, c(10, 0, 0), tolerance = 3e-12)
  expect_equal(x$C_A, c(0, 0, 0), tolerance = 3e-12)
})

test_that("revaluation matches value and derivatives without a positivity promise", {
  for (rate in c(-0.02, 0.05)) {
    maturities <- c(2, 6, 10)
    L <- c(300, 700)
    times <- c(4, 8)
    prices <- 100 / (1 + rate)^maturities
    out <- immunize_duration_convexity(L, times, prices, maturities,
                    maturities * (maturities + 1) / (1 + rate)^2, rate)
    gap <- function(r) sum(out$w * 100 / (1 + r)^maturities) - sum(L / (1 + r)^times)
    h <- 1e-4
    scale <- sum(L / (1 + rate)^times)
    expect_equal(gap(rate) / scale, 0, tolerance = 3e-13)
    expect_lt(abs((gap(rate + h) - gap(rate - h)) / (2 * h * scale)), 1e-5)
    expect_lt(abs((gap(rate + h) - 2 * gap(rate) + gap(rate - h)) / (h^2 * scale)), 1e-4)
    assets <- lapply(maturities, function(t) list(cf = 100, t = t))
    p <- plot_immunization_gap(L, times, assets, out$w, rate, delta = 0.01, n_grid = 21)
    expect_equal(p$data$gap, vapply(p$data$i, gap, numeric(1)), tolerance = 1e-8)
    expect_true(any(p$data$gap > 0) && any(p$data$gap < 0))
  }
})

test_that("gap plots retain weighted cash flows and handle finite cancellations", {
  p <- plot_immunization_gap(100, 1e20, list(list(cf = 110, t = 1e20)),
                            1, 1e-20, delta = 1e-21, n_grid = 5)
  expect_equal(p$data$gap, 10 * exp(-1e20 * log1p(p$data$i)), tolerance = 3e-13)
  p <- plot_immunization_gap(0, 0, list(list(cf = 1e308, t = 1000)),
                            0, -0.5, n_grid = 3)
  expect_equal(p$data$gap, rep(0, 3))
  # Each individual PV overflows, but their matching surplus is exactly zero.
  p <- plot_immunization_gap(1e308, 1000, list(list(cf = 1e308, t = 1000)),
                            1, -0.5, n_grid = 3)
  expect_equal(p$data$gap, rep(0, 3))
  # Weighting must precede discounting when an asset's own PV loses range.
  p <- plot_immunization_gap(0, 0, list(list(cf = 1e308, t = 1000)),
                            1e-310, -0.5, delta = 0.001, n_grid = 3)
  expected <- exp(log(1e308) + log(1e-310) - 1000 * log1p(p$data$i))
  expect_equal(p$data$gap / expected, rep(1, 3), tolerance = 3e-12)
  expect_error(plot_immunization_gap(1e308, 1000, list(list(cf = 0, t = 0)),
                                     1, -0.5, n_grid = 3), "non-finite")
})

test_that("gap grid respects the full annual effective domain and index limits", {
  p <- plot_immunization_gap(100, 1, list(list(cf = 100, t = 1)), 1,
                            -0.9999999, delta = 1e-8, n_grid = 3)
  expect_equal(p$data$gap, rep(0, 3))
  expect_equal(p$data$i, seq(-0.99999991, -0.99999989, length.out = 3), tolerance = 1e-15)
  expect_no_warning(expect_error(plot_immunization_gap(100, 1,
    list(list(cf = 100, t = 1)), 1, 0.05, n_grid = 2^31), "n_grid"))
})
