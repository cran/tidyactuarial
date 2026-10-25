test_that("long bond yields agree with independent discounted payments", {
  for (k in c(1L, 2L, 12L)) {
    periods <- seq_len(30L * k)
    j <- 0.06 / k
    price <- sum((1000 * 0.05 / k) * exp(-periods * log1p(j))) +
      1020 * exp(-30 * k * log1p(j))
    out <- bond_ytm(price, 1000, 0.05, 30, k, R = 1020)
    expect_equal(out$i_period, j, tolerance = 2e-12)
    expect_equal(out$i_effective_annual, expm1(k * log1p(j)),
                 tolerance = 3e-11)
  }
})

test_that("zero-coupon yields agree with the analytic root across currency scales", {
  for (scale in c(1e-250, 1, 1e250)) {
    for (j in c(-0.08, 0, 0.035)) {
      price <- 80 * scale
      redemption <- price * (1 + j)^20
      out <- bond_ytm(price, 100 * scale, 0, 10, k = 2, R = redemption)
      expected <- expm1(log(redemption / price) / 20)
      expect_equal(out$i_period, expected, tolerance = 2e-12)
      expect_equal(out$price, price)
      expect_named(out, c("price", "i_period", "j_nominal",
                          "i_effective_annual", "k"))
    }
  }
  # The automatically chosen domain reaches the representable boundary at -1.
  expect_equal(bond_ytm(1, 1e-8, 0, 1)$i_period, -1 + 1e-8,
               tolerance = 2e-12)
  # A positive cash flow need not have a representable ratio to the price.
  expect_equal(bond_ytm(1e-250, 1e250, 0, 100)$i_period,
               1e5 - 1, tolerance = 1e-10)
  expect_equal(bond_ytm(1e300, 1e-23, 0, 100)$i_period,
               expm1((log(1e-23) - log(1e300)) / 100), tolerance = 2e-12)
})

test_that("bond yield failures do not return unconverged approximations", {
  expect_error(bond_ytm(100, 100, 0.05, 8, interval = c(-0.1, 0.2),
                        maxiter = 1), "converg")
  expect_error(bond_ytm(100, 100, 0.05, 8, interval = c(0.1, Inf)), "interval")
  expect_error(bond_ytm(100, 100, 0.05, 8, tol = Inf), "tol")
  expect_error(bond_ytm(100, 100, 0.05, 8, maxiter = 2^31), "maxiter")
  expect_error(bond_ytm(100, 100, 0.05, 1e-12), "future cash flow")
  expect_error(bond_ytm(100, 100, 0.05, 8, interval = c(0.1, 0.2)),
               "opposite sign|bracket")
  expect_error(bond_ytm(1e-250, 1e250, 0, 1), "representable yield domain")
  expect_error(bond_ytm(1e250, 1e-250, 0, 1), "representable yield domain")
})

test_that("bond prices preserve log yields and representable discounted values", {
  expect_equal(bond_price(1e300, 0, 2, y_effective_per_period = 1e200) /
                 1e-100, 1, tolerance = 2e-13)
  expect_equal(bond_price(1e-300, 0, 100, y_effective_per_period = -0.999) /
                 1, 1, tolerance = 2e-12)
  expect_equal(bond_price(100, 0, 0.5, k = 2, y = -100, y_type = "force") /
                 (100 * exp(50)), 1, tolerance = 2e-13)
  expect_equal(bond_price(1e300, 0, 0.5, k = 2, y = 1000, y_type = "force") /
                 (1e300 * exp(-500)), 1, tolerance = 2e-13)
  expect_equal(bond_price(100, 0, 1e4, k = 1, y = 1e-16) /
                 (100 * exp(-1e-12)), 1, tolerance = 2e-15)
  expect_equal(bond_price(0, 0, 100, y_effective_per_period = -0.999), 0)
  expect_equal(bond_price(100, 0, 1, y_effective_per_period = 0, R = 0), 0)
  expect_error(bond_price(100, 0.05, 8, y = 0.05, tol = Inf), "tol")
  expect_error(bond_price(100, 0.05, 8, y = 0.05, y_m = 2^31), "y_m")
})

test_that("coupon bond prices obey independent sums under all yield families", {
  for (k in c(1L, 2L, 4L)) {
    for (annual in c(-0.1, 0, 0.07)) {
      specs <- list(
        list(y = annual, y_type = "effective"),
        list(y = 4 * ((1 + annual)^0.25 - 1), y_type = "nominal_interest", y_m = 4),
        list(y = 4 * (1 - (1 + annual)^(-0.25)), y_type = "nominal_discount", y_m = 4),
        list(y = log(1 + annual), y_type = "force"))
      times <- seq_len(7 * k) / k
      expected <- sum((500 * 0.06 / k) / (1 + annual)^times) +
        510 / (1 + annual)^7
      for (spec in specs) {
        actual <- do.call(bond_price, c(list(face = 500, c = 0.06,
                                             n = 7, k = k, R = 510), spec))
        expect_equal(actual, expected, tolerance = 2e-13)
      }
    }
  }
})

test_that("bond coupon schedules reject invalid indexing and preserve finite coupons", {
  expect_error(bond_cash_flows(100, 0.05, 1, tol = Inf), "tol")
  expect_error(bond_cash_flows(100, 0.05, 1, k = 2^31), "k")
  expect_error(bond_cash_flows(100, 0.05, 2^31), "n \\* k")
  expect_error(bond_cash_flows(100, 0.05, 1e308, k = 2), "n \\* k")
  expect_error(bond_cash_flows(1e308, 10, 1), "Coupon")
  x <- bond_cash_flows(1e308, 10, 0.1, k = 10, R = 0)
  expect_equal(x$cf, c(1e308, 0), tolerance = 2e-13)
  expect_equal(x$period, c(1L, 1L))
})

test_that("accrued coupon products retain range at coupon and settlement dates", {
  args <- list(face = 1e308, c = 10, k = 10,
               previous_coupon = as.Date("2026-01-01"),
               next_coupon = as.Date("2026-02-01"))
  at_coupon <- do.call(accrued_interest, c(args,
    list(settlement = as.Date("2026-01-01"))))
  expect_equal(at_coupon$coupon_per_period / 1e308, 1, tolerance = 2e-13)
  expect_equal(at_coupon$accrued_interest, 0)
  settled <- do.call(accrued_interest, c(args,
    list(settlement = as.Date("2026-01-16"))))
  expect_equal(settled$accrued_interest / 1e308, 15 / 31, tolerance = 2e-13)
})
