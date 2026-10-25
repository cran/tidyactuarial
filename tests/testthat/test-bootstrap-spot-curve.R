test_that("bootstrap_spot_curve recovers spot discount factors and rates", {
  t <- 1:3
  i <- c(0.04, 0.045, 0.05)
  v <- (1 + i)^(-t)
  bonds <- tibble::tibble(
    instrument = c("Z1", "B2", "B3"),
    price = c(100 * v[1],
              5 * v[1] + 105 * v[2],
              6 * v[1] + 6 * v[2] + 106 * v[3]),
    t = list(1, c(1, 2), c(1, 2, 3)),
    cf = list(100, c(5, 105), c(6, 6, 106))
  )

  out <- bootstrap_spot_curve(bonds)
  expect_s3_class(out, "tbl_df")
  expect_equal(out$instrument, c("Z1", "B2", "B3"))
  expect_equal(out$t, as.numeric(t))
  expect_equal(out$i, i, tolerance = 1e-12)
  expect_equal(out$v, v, tolerance = 1e-12)
  expect_equal(out$repriced_price, bonds$price, tolerance = 1e-12)
  expect_lt(max(abs(out$repricing_error)), 1e-10)

  existing <- yield_curve(t = out$t, i = out$i)
  expect_equal(existing$v[[1]], out$v, tolerance = 1e-12)
  expect_equal(
    forward_rate(t = out$t, i = out$i, t_start = 1, t_end = 2)$f,
    (v[1] / v[2]) - 1,
    tolerance = 1e-12
  )
})

test_that("bootstrap_spot_curve sorts instruments and supports custom columns", {
  v <- (1 + c(0.04, 0.05))^(-(1:2))
  x <- tibble::tibble(
    name = c("B2", "Z1"),
    market = c(5 * v[1] + 105 * v[2], 100 * v[1]),
    dates = list(c(2, 1), 1),
    payments = list(c(105, 5), 100)
  )
  out <- bootstrap_spot_curve(
    x, col_instrument = "name", col_price = "market",
    col_t = "dates", col_cf = "payments"
  )
  expect_equal(out$t, c(1, 2))
  expect_equal(out$i, c(0.04, 0.05), tolerance = 1e-12)
})

test_that("bootstrap_spot_curve combines simultaneous coupon and redemption", {
  v <- (1 + c(0.04, 0.05))^(-(1:2))
  x <- tibble::tibble(
    instrument = c("Z1", "B2"),
    price = c(100 * v[1], 5 * v[1] + 105 * v[2]),
    t = list(1, c(1, 2, 2)),
    cf = list(100, c(5, 5, 100))
  )
  out <- bootstrap_spot_curve(x)
  expect_equal(out$terminal_cf, c(100, 105))
  expect_equal(out$i, c(0.04, 0.05), tolerance = 1e-12)
})

test_that("bootstrap_spot_curve handles zero coupons and negative spot rates", {
  v <- (1 + c(-0.01, -0.02))^(-(1:2))
  x <- tibble::tibble(
    instrument = c("Z1", "Z2"),
    price = 100 * v,
    t = list(1, c(1, 2)),
    cf = list(100, c(0, 100))
  )
  out <- bootstrap_spot_curve(x)
  expect_equal(out$i, c(-0.01, -0.02), tolerance = 1e-12)
  expect_true(all(out$v > 1))
})

test_that("bootstrap_spot_curve rejects an unidentified earlier node", {
  x <- tibble::tibble(
    instrument = c("Z1", "B3"),
    price = c(96, 100),
    t = list(1, c(1, 2, 3)),
    cf = list(100, c(5, 5, 105))
  )
  expect_error(bootstrap_spot_curve(x), "without a bootstrapped discount factor")
})

test_that("bootstrap_spot_curve rejects invalid implied discounts", {
  x <- tibble::tibble(
    instrument = c("Z1", "B2"),
    price = c(95, 1),
    t = list(1, c(1, 2)),
    cf = list(100, c(10, 100))
  )
  expect_error(bootstrap_spot_curve(x), "nonpositive")
})

test_that("bootstrap_spot_curve rejects duplicate final maturities", {
  x <- tibble::tibble(
    instrument = c("A", "B"),
    price = c(95, 97),
    t = list(1, 1),
    cf = list(100, 100)
  )
  expect_error(bootstrap_spot_curve(x), "distinct final maturity")
})

test_that("bootstrap_spot_curve validates observed inputs", {
  x <- tibble::tibble(instrument = "A", price = 95, t = list(1), cf = list(100))
  expect_error(bootstrap_spot_curve(x[0, ]), "nonempty data frame")
  expect_error(bootstrap_spot_curve(x, col_cf = "missing"), "Missing required")
  expect_error(bootstrap_spot_curve(dplyr::bind_rows(x, x)), "unique")
  bad_price <- x
  bad_price$price <- NA_real_
  expect_error(bootstrap_spot_curve(bad_price), "finite positive")
  bad_cf <- x
  bad_cf$cf <- list(-100)
  expect_error(bootstrap_spot_curve(bad_cf), "Invalid payment")
  bad_t <- x
  bad_t$t <- list(0)
  expect_error(bootstrap_spot_curve(bad_t), "Invalid payment")
  bad_cf2 <- x
  bad_cf2$cf <- list(c(0, 0))
  bad_cf2$t <- list(c(1, 2))
  expect_error(bootstrap_spot_curve(bad_cf2), "positive cash flow")
})

test_that("bootstrap_spot_curve accepts bond_cash_flows schedules directly", {
  z <- bond_cash_flows(face = 100, c = 0, n = 1, k = 1)
  b <- bond_cash_flows(face = 100, c = 0.05, n = 2, k = 1)
  i <- c(0.04, 0.05)
  v <- (1 + i)^(-(1:2))
  x <- tibble::tibble(
    instrument = c("zero", "coupon"),
    price = c(100 * v[1], 5 * v[1] + 105 * v[2]),
    t = list(z$t, b$t),
    cf = list(z$cf, b$cf)
  )
  expect_equal(bootstrap_spot_curve(x)$i, i, tolerance = 1e-12)
})
