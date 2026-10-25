test_that("mc_annuity makes per-payment and annualized units explicit", {
  simulations <- tibble::tibble(
    Kx = c(2, 2),
    Tx = c(2.4, 2.4)
  )

  result <- mc_annuity(
    simulations,
    i = 0.05,
    payment = 1 / 12,
    k = 12,
    type = "whole",
    timing = "due"
  )

  expect_equal(result$payment_per_payment, rep(1 / 12, 2))
  expect_equal(result$payment_annualized, rep(1, 2))
  expect_equal(result$payment, result$payment_per_payment)
})


test_that("mc_annuity annual whole-life due payment times are correct", {
  simulations <- tibble::tibble(
    Kx = c(0, 2, 4)
  )

  result <- mc_annuity(
    simulations,
    i = 0,
    payment = 10,
    k = 1,
    type = "whole",
    timing = "due"
  )

  expect_equal(result$n_payments, c(1L, 3L, 5L))
  expect_equal(result$pv_annuity, c(10, 30, 50))
})


test_that("mc_annuity fractional temporary schedule uses payment of 1 over k", {
  simulations <- tibble::tibble(
    Kx = 3,
    Tx = 3.2
  )

  result <- mc_annuity(
    simulations,
    i = 0,
    payment = 1 / 4,
    k = 4,
    type = "temporary",
    n = 2,
    timing = "due"
  )

  expect_equal(result$n_payments, 8L)
  expect_equal(result$pv_annuity, 2)
  expect_equal(result$last_payment_time, 1.75)
})


test_that("mc_annuity requires finite terms to align with the payment grid", {
  simulations <- tibble::tibble(
    Kx = 3,
    Tx = 3.2
  )

  expect_error(
    mc_annuity(
      simulations,
      i = 0.05,
      k = 4,
      type = "temporary",
      n = 2.3,
      timing = "due"
    ),
    "`n \\* k` must be an integer"
  )

  expect_error(
    mc_annuity(
      simulations,
      i = 0.05,
      k = 12,
      type = "guaranteed",
      n_guar = 2.3,
      timing = "due"
    ),
    "`n_guar \\* k` must be an integer"
  )
})


test_that("mc_annuity rejects noninteger interest conversion frequency", {
  simulations <- tibble::tibble(
    Kx = 2
  )

  expect_error(
    mc_annuity(
      simulations,
      i = 0.06,
      i_type = "nominal_interest",
      m = 1.5,
      type = "whole"
    ),
    "`m`"
  )
})


test_that("mc_annuity certain payments do not require complete lifetimes", {
  simulations <- tibble::tibble(
    Kx = c(NA_real_, NA_real_)
  )

  result <- mc_annuity(
    simulations,
    i = 0,
    payment = 100,
    k = 12,
    type = "certain",
    n = 1,
    timing = "immediate"
  )

  expect_equal(result$n_payments, c(12L, 12L))
  expect_equal(result$pv_annuity, c(1200, 1200))
})


test_that("mc_annuity deprecated aliases use exact conflict detection", {
  simulations <- tibble::tibble(
    Kx = 2
  )

  expect_error(
    mc_annuity(
      simulations,
      i = 0.05,
      k = 1,
      payments_per_year = 1,
      type = "whole"
    ),
    "Provide only one of `k`"
  )

  result <- mc_annuity(
    data = simulations,
    rate = 0.05,
    payments_per_year = 1,
    annuity = "whole_life"
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(result$type, "whole")
})
test_that("missing lifetimes do not become zero-value life annuities", {
  annual <- mc_annuity(tibble::tibble(Kx = c(2, NA_real_)), i = 0)
  expect_true(is.na(annual$pv_annuity[2]))
  expect_true(is.na(annual$n_payments[2]))
  fractional <- mc_annuity(tibble::tibble(Kx = c(2, 2), Tx = c(2.4, NA_real_)),
                           i = 0, k = 4)
  expect_true(is.na(fractional$pv_annuity[2]))
  expect_true(is.na(fractional$n_payments[2]))
})

test_that("fractional payments distinguish deaths on either side of a payment date", {
  lifetimes <- c(1e-10, .25 - 1e-10, .25, .25 + 1e-10)
  data <- tibble::tibble(Kx = floor(lifetimes), Tx = lifetimes)
  due <- mc_annuity(data, i = 0, payment = 1, k = 4, timing = "due")
  immediate <- mc_annuity(data, i = 0, payment = 1, k = 4, timing = "immediate")
  # Preserve the existing endpoint conventions: due dates < T;
  # immediate dates <= T. Nearby distinct lifetimes must remain distinct.
  expect_equal(due$n_payments, c(1L, 1L, 1L, 2L))
  expect_equal(immediate$n_payments, c(0L, 0L, 1L, 1L))
  reserve <- mc_reserve(data, t = 0, i = 0, P = 4, premium_unit = "annualized",
    payment = 1, k = 4, premium_timing = "due", benefit = 100)
  expect_equal(reserve$Y_t, due$pv_annuity)
  expect_equal(reserve$L_t, 100 - due$pv_annuity)
})

test_that("guaranteed fractional annuities pay through the guarantee after immediate death", {
  data <- tibble::tibble(Kx = 0, Tx = 0)
  for (timing in c("due", "immediate")) {
    out <- mc_annuity(data, i = 0, payment = 1, k = 4,
      type = "guaranteed", n_guar = 1, timing = timing)
    expect_equal(out$n_payments, 4L)
    expect_equal(out$pv_annuity, 4)
  }
})
test_that("fractional deferred calendars are anchored at the deferment", {
  T <- c(.1, .2, .375 - 1e-10, .375, .375 + 1e-10, .7)
  data <- data.frame(Kx = floor(T), Tx = T)
  for (type in c("deferred", "deferred_temporary")) {
    for (timing in c("due", "immediate")) {
      # Explicit reference grid, shifted by h=1/8 rather than issue time.
      grid <- if (timing == "due") c(.125, .375, .625) else c(.375, .625)
      if (type == "deferred_temporary" && timing == "due") grid <- head(grid, 2)
      expected <- vapply(T, function(t) {
        dates <- if (timing == "due") grid[grid < t] else grid[grid <= t]
        sum(1.05^(-dates))
      }, numeric(1))
      out <- mc_annuity(data, i = .05, payment = 1, k = 4, h = .125,
        type = type, n = if (type == "deferred_temporary") .5 else NULL,
        timing = timing)
      expect_equal(out$pv_annuity, expected, tolerance = 1e-12)
      reserve <- mc_reserve(data, t = 0, i = .05, P = 4,
        premium_unit = "annualized", payment = 1, k = 4, h = .125,
        annuity_type = type, n = if (type == "deferred_temporary") .5 else NULL,
        premium_timing = timing)
      expect_equal(reserve$Y_t, expected, tolerance = 1e-12)
    }
  }
})
