test_that("forward rates retain flat and nonflat log growth without powers", {
  expect_equal(forward_rate(t = c(10000, 20000), i = c(0.1, 0.1),
                             t_start = 10000, t_end = 20000)$f, 0.1,
               tolerance = 2e-14)
  expect_equal(forward_rate(t = c(1, 2), i = c(1e-20, 1e-20),
                             t_start = 1, t_end = 2)$f / 1e-20, 1,
               tolerance = 2e-14)
  expect_equal(forward_rate(t = c(0.01, 0.02), i = c(-100, -100),
                             i_type = "force", t_start = 0.01, t_end = 0.02)$f,
               expm1(-100))
  expected <- expm1(0.06 + (0.06 - 0.04) / 2)
  out <- forward_rate(t = c(1, 3), i = c(0.04, 0.06), i_type = "force",
                        t_start = 1, t_end = 3)
  expect_equal(out$f, expected, tolerance = 2e-14)
  expect_equal(out$i_start, expm1(0.04), tolerance = 2e-14)
  # A large common origin must not erase a flat interval growth rate.
  expect_equal(forward_rate(t = c(1e12, 1e12 + 1), i = c(0.03, 0.03),
                             i_type = "force", t_start = 1e12,
                             t_end = 1e12 + 1)$f, expm1(0.03), tolerance = 2e-14)
})

test_that("linear forward extraction still interpolates annual effective rates", {
  times <- c(10000, 20000, 30000)
  rates <- c(0.03, 0.05, 0.06)
  out <- forward_rate(t = times, i = rates, t_start = 15000, t_end = 25000,
                        method = "linear")
  expect_equal(out$i_start, 0.04, tolerance = 2e-14)
  expect_equal(out$i_end, 0.055, tolerance = 2e-14)
  expected <- expm1((25000 * log(1.055) - 15000 * log(1.04)) / 10000)
  expect_equal(out$f, expected, tolerance = 2e-14)
})

test_that("adjacent forwards preserve small original log rates and composition", {
  expect_equal(forward_curve(t = c(0, 1, 2), i = rep(1e-20, 3))$f / 1e-20,
               c(1, 1), tolerance = 2e-14)
  out <- forward_curve(t = c(3, 1, 2), i = c(0.06, 0.04, 0.05),
                        i_type = "force", include_origin = TRUE)
  expect_equal(out$f, expm1(c(0.04, 0.06, 0.08)), tolerance = 2e-14)
  expect_equal(cumprod((1 + out$f)^out$duration), exp(c(0.04, 0.1, 0.18)),
               tolerance = 2e-14)
  out <- forward_curve(t = c(1e12, 1e12 + 1), i = c(1e-12, 1e-12),
                        i_type = "force")
  expect_equal(out$f / expm1(1e-12), c(1, 1), tolerance = 2e-14)
})

test_that("yield discounts retain the original basis and zero-time identity", {
  out <- yield_curve(t = c(0, 0.01, 0.02), i = rep(-100, 3), i_type = "force")
  expect_equal(out$v[[1]], exp(c(0, 1, 2)), tolerance = 2e-14)
  expect_equal(out$i_effective[[1]], rep(-1, 3))
  expect_equal(yield_curve(t = c(1e12, 2e12), i = c(1e-20, 1e-20))$v[[1]],
               exp(-c(1e-8, 2e-8)), tolerance = 2e-15)
})

test_that("missing curves keep plot slots and endpoint propagation", {
  for (missing in list(c(1L), c(2L), c(1L, 2L))) {
    d <- tibble::tibble(id = c("a", "b"), t = list(1:2, 1:2),
                       i = list(c(0.04, 0.05), c(0.04, 0.05)))
    for (a in missing) d$i[[a]][1] <- NA_real_
    out <- yield_curve(d, plot = TRUE)
    expect_equal(length(out$yield_curve_plot), 2L)
    for (a in 1:2) expect_identical(is.null(out$yield_curve_plot[[a]]), a %in% missing)
    d$t_start <- c(1, 1)
    d$t_end <- c(2, 2)
    out <- forward_rate(d, plot = TRUE)
    expect_identical(is.na(out$f), 1:2 %in% missing)
    for (a in 1:2) expect_identical(is.null(out$forward_rate_plot[[a]]), a %in% missing)
  }
  d <- tibble::tibble(t = list(1:2), i = list(c(0.04, 0.05)),
                     t_start = NA_real_, t_end = 2)
  expect_true(is.na(forward_rate(d)$f))
  expect_error(forward_rate(d, .na = "error"), "Missing values")
  expect_warning(dropped <- forward_rate(d, .na = "drop", plot = TRUE), "Dropping")
  expect_equal(nrow(dropped), 0L)
  expect_equal(length(dropped$forward_rate_plot), 0L)
})

test_that("par coupons are independent of currency scale and retain small rates", {
  for (face in c(1e-250, 100, 1e308)) {
    out <- par_rate(t = 1:2, i = c(0.05, 0.05), face = face)
    expect_equal(out$c_par, 0.05, tolerance = 3e-13)
    expect_equal(out$coupon_per_period / face, 0.05, tolerance = 3e-13)
    expect_equal(out$repriced_price / face, 1, tolerance = 3e-13)
  }
  for (rate in c(-1e-20, 0, 1e-20)) {
    out <- par_rate(t = 1:2, i = rep(rate, 2))
    if (rate == 0) expect_equal(out$c_par, 0) else
      expect_equal(out$c_par / rate, 1, tolerance = 3e-13)
  }
  out <- par_rate(t = c(0.5, 1), i = c(0.04, 0.04), k = 2,
                    face = 1e308, R = 0.9e308)
  v <- 1.04^(-c(0.5, 1))
  expected <- 2 * (1 - 0.9 * v[2]) / sum(v)
  expect_equal(out$c_par, expected, tolerance = 3e-13)
  expect_equal(out$repriced_price / 1e308, 1, tolerance = 3e-13)
})

test_that("bootstrap recovers a triangular system under currency changes and ordering", {
  for (scale in c(1e-200, 1, 1e200)) {
    discounts <- c(1.02, 0.95, 0.87)
    payments <- rbind(c(100, 0, 0), c(4, 104, 0), c(6, 6, 106)) * scale
    data <- tibble::tibble(instrument = c("a", "b", "c"),
      price = as.vector(payments %*% discounts),
      t = list(1, 1:2, c(3, 1, 2, 3)),
      cf = list(payments[1, 1], payments[2, 1:2], c(100, 6, 6, 6) * scale))
    out <- bootstrap_spot_curve(data[c(3, 1, 2), ])
    expect_equal(out$instrument, c("a", "b", "c"))
    expect_equal(out$v, discounts, tolerance = 3e-13)
    expect_equal(out$i, discounts^(-1 / (1:3)) - 1, tolerance = 3e-13)
    expect_equal(out$repriced_price / out$price, rep(1, 3), tolerance = 3e-13)
    par <- out |> par_rate(face = 100)
    expect_equal(par$c_par, (1 - discounts[3]) / sum(discounts), tolerance = 3e-13)
    expect_equal(par$coupon_per_period * sum(discounts) + 100 * discounts[3],
                 100, tolerance = 3e-13)
    forwards <- out |> forward_curve(include_origin = TRUE)
    expect_equal(cumprod((1 + forwards$f)^forwards$duration), 1 / discounts,
                 tolerance = 3e-13)
  }
})

test_that("nodewise conventions follow the same sorted or list-column curve", {
  times <- c(1, 2, 4)
  growth <- c(0.02, 0.03, 0.04)
  rates <- c(2 * expm1(0.02 / 2), 0.03, -4 * expm1(-0.04 / 4))
  types <- c("nominal_interest", "force", "nominal_discount")
  frequency <- c(2, 1, 4)
  y <- yield_curve(t = times, i = rates, i_type = types, m = frequency)
  expect_equal(y$v[[1]], exp(-times * growth), tolerance = 3e-14)
  f <- forward_rate(t = times, i = rates, i_type = types, m = frequency,
                      t_start = 2, t_end = 4)
  expect_equal(f$f, expm1(0.05), tolerance = 3e-14)
  order <- c(3, 1, 2)
  f <- forward_curve(t = times[order], i = rates[order],
                       i_type = types[order], m = frequency[order], include_origin = TRUE)
  expect_equal(f$f, expm1(c(0.02, 0.04, 0.05)), tolerance = 3e-14)
  curves <- tibble::tibble(id = c("first", "second"), t = list(1:3, 1:3),
                          i = list(rep(0.04, 3), rep(0.06, 3)))
  result <- curves |> par_rate(face = c(100, 200))
  expect_identical(result$id, curves$id)
  expect_equal(result$c_par, c(0.04, 0.06), tolerance = 3e-13)
})
