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
