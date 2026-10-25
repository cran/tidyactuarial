test_that("book values obey prospective sums, retrospective values and the ledger", {
  for (k in c(1L, 2L, 4L)) {
    for (j in c(-0.05, 0, 0.04)) {
      N <- 3L * k
      C <- 500 * 0.06 / k
      pay <- rep(C, N)
      pay[N] <- pay[N] + 510
      periods <- seq_len(N)
      price <- sum(pay / (1 + j)^periods)
      expected <- vapply(0:N, function(s) {
        future <- periods > s
        sum(pay[future] / (1 + j)^(periods[future] - s))
      }, numeric(1))
      actual <- bond_book_value(500, 0.06, 3, t = (0:N) / k, k = k,
                                 y_effective_per_period = j, R = 510)
      expect_equal(actual, expected, tolerance = 2e-13)
      expect_equal(actual[-1], actual[-length(actual)] * (1 + j) - pay,
                   tolerance = 2e-12)
      retro <- vapply(0:N, function(s) {
        past <- periods <= s
        price * (1 + j)^s - sum(pay[past] * (1 + j)^(s - periods[past]))
      }, numeric(1))
      expect_equal(actual, retro, tolerance = 2e-12)
      shuffled <- c(N, 0L, 1L, 1L)
      out <- bond_book_value(500, 0.06, 3, t = shuffled / k, k = k,
                              y_effective_per_period = j, R = 510, tidy = TRUE)
      expect_equal(out$book_value, expected[shuffled + 1L], tolerance = 2e-13)
      expect_equal(out$valuation_period, shuffled)
      expect_named(out, c("t", "valuation_period", "book_value", "yield_per_period",
                          "yield_effective_annual", "k", "face", "c", "n", "R"))
    }
  }
})

test_that("book values retain original annual log yields and monetary range", {
  out <- bond_book_value(100, 0, 1, t = c(0, 0.5, 1), k = 2,
                          y = -100, y_type = "force", tidy = TRUE)
  expect_equal(out$book_value[1:2] / (100 * exp(c(100, 50))), c(1, 1),
               tolerance = 2e-13)
  expect_equal(out$book_value[3], 0)
  expect_equal(out$yield_effective_annual, rep(-1, 3))
  expect_equal(bond_book_value(1e300, 0, 2, t = c(0, 1, 2),
                                y_effective_per_period = 1e200) /
                 c(1e-100, 1e100, 1), c(1, 1, 0), tolerance = 2e-13)
  expect_equal(bond_book_value(1e300, 0, 0.5, t = 0, k = 2,
                                y = 1000, y_type = "force") /
                 (1e300 * exp(-500)), 1, tolerance = 2e-13)
  expect_equal(bond_book_value(1e300, 0, 1, t = 0,
                                y = 1000, y_type = "force") /
                 exp(log(1e300) - 1000), 1, tolerance = 2e-13)
  expect_equal(bond_book_value(100, 0, 1e4, t = 0, y = 1e-16) /
                 (100 * exp(-1e-12)), 1, tolerance = 2e-15)
  expect_equal(bond_book_value(100, 0, 0, t = c(0, 0), y = 0.05), c(0, 0))
  empty <- bond_book_value(100, 0.05, 2, t = numeric(), y = 0.05, tidy = TRUE)
  expect_equal(nrow(empty), 0L)
  expect_type(empty$valuation_period, "integer")
  expect_equal(bond_book_value(0, 0, 2, t = 0:2,
                                y_effective_per_period = -0.999), rep(0, 3))
})

test_that("call scenarios pair redemptions with dates and include the final coupon", {
  for (k in c(1L, 2L, 4L)) {
    for (j in c(-0.03, 0, 0.05)) {
      times <- c(2, 1, 3)
      redemptions <- c(515, 520, 510)
      expected <- vapply(seq_along(times), function(a) {
        periods <- seq_len(times[a] * k)
        sum((500 * 0.06 / k) / (1 + j)^periods) +
          redemptions[a] / (1 + j)^(times[a] * k)
      }, numeric(1))
      out <- bond_callable_price(500, 0.06, 3, k, call_t = times[1:2],
                                  call_R = redemptions[1:2], R = 510,
                                  y_effective_per_period = j, tidy = TRUE)
      expect_equal(out$scenario_id, c(2L, 1L, 3L))
      expect_equal(out$redemption_value, c(520, 515, 510))
      expect_equal(out$price_at_target_yield, expected[c(2, 1, 3)], tolerance = 2e-13)
      expect_equal(bond_callable_price(500, 0.06, 3, k, call_t = times[1:2],
                                        call_R = redemptions[1:2], R = 510,
                                        y_effective_per_period = j), min(expected),
                   tolerance = 2e-13)
      expect_named(out, c("scenario_id", "scenario_type", "t", "redemption_value",
                          "price_at_target_yield", "is_worst_case", "yield_per_period",
                          "yield_effective_annual", "k", "face", "c", "n", "R"))
    }
  }
})

test_that("call scenarios preserve the original yield specification", {
  for (force in c(-100, 1000)) {
    face <- if (force > 0) 1e300 else 100
    out <- bond_callable_price(face, 0, 1, k = 2, call_t = 0.5, call_R = face,
                                 y = force, y_type = "force", tidy = TRUE)
    expected <- face * exp(-force * c(0.5, 1))
    # The final scenario at force 1000 is recovered even though exp(-1000) is zero.
    if (force > 0) expected[2] <- exp(log(face) - force)
    expect_equal(out$price_at_target_yield / expected, c(1, 1), tolerance = 3e-13)
  }
  expect_equal(bond_callable_price(100, 0.05, 2, call_t = 1, call_R = 100,
                                    y_effective_per_period = 0.04, y = -2),
                 105 / 1.04, tolerance = 2e-13)
  saturated <- bond_callable_price(1e300, 0, 2, call_t = 1, call_R = 1e300,
                                     y = 1000, y_type = "force", tidy = TRUE)
  expect_equal(saturated$price_at_target_yield[1] / exp(log(1e300) - 1000),
               1, tolerance = 2e-13)
  expect_equal(saturated$price_at_target_yield[2], 0)
  expect_equal(saturated$yield_per_period, rep(Inf, 2))
})

test_that("annual bases preserve the same remaining payments and call scenarios", {
  annual <- 0.07
  specs <- list(list(y = annual, y_type = "effective"),
                list(y = 4 * ((1 + annual)^0.25 - 1), y_type = "nominal_interest", y_m = 4),
                list(y = 4 * (1 - (1 + annual)^(-0.25)), y_type = "nominal_discount", y_m = 4),
                list(y = log(1 + annual), y_type = "force"))
  for (spec in specs) {
    book <- do.call(bond_book_value,
      c(list(face = 100, c = 0.05, n = 2, t = c(0, 1, 2), k = 2), spec))
    expected <- c(sum(2.5 / 1.07^((1:4) / 2)) + 100 / 1.07^2,
                  2.5 / 1.07^0.5 + 102.5 / 1.07, 0)
    expect_equal(book, expected, tolerance = 2e-13)
    callable <- do.call(bond_callable_price,
      c(list(face = 100, c = 0.05, n = 2, k = 2, call_t = 1, call_R = 102), spec))
    expect_equal(callable, min(2.5 / 1.07^0.5 + 104.5 / 1.07, expected[1]),
                 tolerance = 2e-13)
  }
  data <- tibble::tibble(id = c("initial", "after_coupon", "paid"), t = 0:2)
  out <- data |> dplyr::mutate(value = bond_book_value(100, 0.05, 2, t, y = 0))
  expect_identical(out$id, data$id)
  expect_equal(out$value, c(110, 105, 0))
})

test_that("C05 identifies computed minimum prices independently of monetary units", {
  for (scale in c(1e-12, 1, 1e12)) {
    out <- bond_callable_price(100 * scale, 0.08, 10, call_t = c(5, 8),
                                 call_R = c(102, 101) * scale, y = 0.05, tidy = TRUE)
    expect_identical(out$is_worst_case, c(TRUE, FALSE, FALSE))
    expect_identical(out$is_worst_case,
                     out$price_at_target_yield == min(out$price_at_target_yield))
  }
  # At zero yield these three computed scenario prices are exactly equal.
  tie <- bond_callable_price(100, 0.02, 3, call_t = c(1, 2),
                               call_R = c(104, 102), R = 100, y = 0, tidy = TRUE)
  expect_equal(tie$price_at_target_yield, rep(106, 3))
  expect_identical(tie$is_worst_case, rep(TRUE, 3))
  # Nearly equal prices are not computed ties; tol is a time-alignment tolerance.
  close <- bond_callable_price(100, 0, 2, call_t = 1, call_R = 100 + 1e-6,
                                 R = 100, y = 0, tidy = TRUE)
  expect_identical(close$is_worst_case, c(FALSE, TRUE))
})

test_that("book and call schedules reject ambiguous or unrepresentable period indices", {
  expect_no_warning(expect_error(bond_book_value(100, 0.05, 2^31, t = 0,
                                                  y = 0.05), "n \\* k"))
  expect_error(bond_book_value(100, 0.05, 1e308, t = 0, k = 2, y = 0.05), "n \\* k")
  expect_error(bond_callable_price(100, 0.05, 1e308, k = 2, call_t = 1,
                                    call_R = 100, y = 0.05), "n \\* k")
  expect_error(bond_callable_price(100, 0.05, 5, call_t = c(1, 1 + 1e-12),
                                    call_R = c(100, 110), y = 0.05), "coupon period")
  expect_error(bond_callable_price(100, 0.05, 5, call_t = 1e-12,
                                    call_R = 100, y = 0.05), "coupon period")
  expect_error(bond_callable_price(100, 0.05, 5, call_t = 5 - 1e-12,
                                    call_R = 100, y = 0.05), "coupon period")
})
