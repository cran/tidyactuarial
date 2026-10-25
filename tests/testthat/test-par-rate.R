test_that("par_rate recovers a nonflat annual par coupon rate and reprices at face", {
  t <- 1:3
  i <- c(0.04, 0.045, 0.05)
  v <- (1 + i)^(-t)
  x <- par_rate(t = t, i = i, face = 100)
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), 1L)
  expect_equal(x$c_par, (1 - v[[3]]) / sum(v), tolerance = 1e-12)
  expect_equal(x$annuity_factor, sum(v), tolerance = 1e-12)
  expect_equal(x$v_n, v[[3]], tolerance = 1e-12)
  expect_equal(x$repriced_price, 100, tolerance = 1e-12)
  expect_equal(x$repricing_error, 0, tolerance = 1e-12)
})

test_that("par_rate accepts the nodewise output of bootstrap_spot_curve", {
  i <- c(0.04, 0.045, 0.05)
  v <- (1 + i)^(-(1:3))
  instruments <- tibble::tibble(
    instrument = c("Z1", "B2", "B3"),
    price = c(100 * v[[1L]], 5 * v[[1L]] + 105 * v[[2L]],
              6 * v[[1L]] + 6 * v[[2L]] + 106 * v[[3L]]),
    t = list(1, c(1, 2), c(1, 2, 3)),
    cf = list(100, c(5, 105), c(6, 6, 106))
  )
  x <- bootstrap_spot_curve(instruments) |> par_rate()
  expect_equal(x$c_par, (1 - v[[3L]]) / sum(v), tolerance = 1e-12)
  expect_equal(x$repriced_price, 100, tolerance = 1e-12)
})

test_that("par_rate supports semiannual coupons and reproduces bond flows", {
  t <- seq_len(4L) / 2
  i <- rep(0.04, length(t))
  x <- par_rate(t = t, i = i, n = 2, face = 1000, k = 2)
  expected <- 2 * (sqrt(1.04) - 1)
  expect_equal(x$c_par, expected, tolerance = 1e-12)
  flows <- bond_cash_flows(face = 1000, c = x$c_par, n = 2, k = 2)
  present_value <- sum(flows$cf * discount_factor_spot(t = flows$t, i = 0.04))
  expect_equal(present_value, 1000, tolerance = 1e-9)
})

test_that("par_rate accounts for redemption different from face value", {
  t <- 1:2
  i <- c(0.04, 0.05)
  v <- (1 + i)^(-t)
  x <- par_rate(t = t, i = i, face = 100, R = 103)
  expect_equal(x$c_par, (100 - 103 * v[[2L]]) /
                 (100 * sum(v)), tolerance = 1e-12)
  expect_equal(x$repriced_price, 100, tolerance = 1e-12)
})

test_that("par_rate accepts yield_curve output and preserves multi-curve IDs", {
  x <- yield_curve(t = 1:3, i = c(0.04, 0.045, 0.05)) |>
    par_rate()
  expect_equal(nrow(x), 1L)
  expect_true(is.numeric(x$c_par))

  curves <- tibble::tibble(
    curve_id = c("A", "B"),
    t = list(1:2, 1:2),
    i = list(c(0.04, 0.05), c(0.03, 0.035))
  )
  y <- curves |> par_rate()
  expect_equal(y$curve_id, c("A", "B"))
  expect_equal(y$c_par,
               c((1 - 1.05^(-2)) / (1.04^(-1) + 1.05^(-2)),
                 (1 - 1.035^(-2)) / (1.03^(-1) + 1.035^(-2))),
               tolerance = 1e-12)
})

test_that("par_rate works with force of interest and custom input columns", {
  x <- tibble::tibble(maturity = 1:2, spot = rep(log1p(0.04), 2))
  y <- par_rate(x, col_t = "maturity", col_i = "spot", i_type = "force")
  z <- par_rate(t = 1:2, i = rep(0.04, 2))
  expect_equal(y$c_par, z$c_par, tolerance = 1e-12)
})

test_that("par_rate requires each coupon date to be known", {
  expect_error(par_rate(t = c(1, 3), i = c(0.04, 0.05), k = 1),
               "no unique discount-factor node")
  expect_error(par_rate(t = 1:2, i = c(0.04, 0.05), k = 2),
               "no unique discount-factor node")
  expect_error(par_rate(t = c(1, 1, 2), i = c(0.04, 0.05, 0.06)),
               "duplicate or indistinguishable")
  expect_error(par_rate(t = 1:2, i = c(0.04, 0.05), n = 1.5),
               "n \\* k")
})

test_that("par_rate validates inputs and handles negative spot rates", {
  expect_error(par_rate(t = c(1, 2), i = c(0.04)), "equal nonzero length")
  expect_error(par_rate(t = 1:2, i = c(0.04, -1)), "greater than -1")
  expect_error(par_rate(t = 1:2, i = c(0.04, NA_real_)), "finite numeric")
  expect_error(par_rate(t = 1:2, i = c(0.04, 0.05), face = 0), "strictly positive")
  expect_error(par_rate(t = 1:2, i = c(0.04, 0.05), k = 1.5), "positive integers")
  expect_error(par_rate(t = 1:2, i = c(0.04, 0.05), R = -1), "nonnegative")
  neg <- par_rate(t = 1:2, i = rep(-0.01, 2))
  expect_lt(neg$c_par, 0)
  expect_equal(neg$repriced_price, 100, tolerance = 1e-12)
})

test_that("par_rate uses maturity and financial parameters rowwise", {
  curves <- tibble::tibble(
    curve_id = c("One", "Two"),
    t = list(1, 1:2),
    i = list(0.04, c(0.03, 0.04))
  )
  x <- par_rate(curves, n = c(1, 2), k = 1, face = c(100, 1000))
  expect_equal(x$n, c(1, 2))
  expect_equal(x$face, c(100, 1000))
  expect_equal(x$c_par[[1]], 0.04, tolerance = 1e-12)
  expect_equal(x$repriced_price, c(100, 1000), tolerance = 1e-10)
})
