test_that("mc_reserve estimates an annualized premium on a normalized monthly basis", {
  simulated <- tibble::tibble(
    Kx = c(1, 2),
    Tx = c(1.5, 2.5)
  )

  result <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    benefit = 100,
    payment = 1 / 12,
    k = 12,
    type = "whole",
    annuity_type = "whole",
    timing = "end_of_year",
    premium_timing = "due",
    reserve_timing = "before_payment"
  )

  expect_equal(
    unique(result$premium_annualized),
    50,
    tolerance = 1e-12
  )

  expect_equal(
    unique(result$premium_per_payment),
    50 / 12,
    tolerance = 1e-12
  )

  expect_equal(
    unique(result$premium_annuity_scale),
    50,
    tolerance = 1e-12
  )

  expect_lt(abs(mean(result$L_t)), 1e-12)
})


test_that("mc_reserve is invariant to an equivalent premium-annuity basis", {
  simulated <- tibble::tibble(
    Kx = c(1, 2),
    Tx = c(1.5, 2.5)
  )

  normalized <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    benefit = 100,
    payment = 1 / 12,
    k = 12,
    type = "whole",
    annuity_type = "whole"
  )

  unit_payments <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    benefit = 100,
    payment = 1,
    k = 12,
    type = "whole",
    annuity_type = "whole"
  )

  expect_equal(
    normalized$premium_annualized,
    unit_payments$premium_annualized,
    tolerance = 1e-12
  )

  expect_equal(
    normalized$premium_per_payment,
    unit_payments$premium_per_payment,
    tolerance = 1e-12
  )

  expect_equal(
    normalized$L_t,
    unit_payments$L_t,
    tolerance = 1e-12
  )

  expect_equal(
    unique(unit_payments$premium_annuity_scale),
    50 / 12,
    tolerance = 1e-12
  )
})


test_that("mc_reserve distinguishes before- and after-payment reserves", {
  simulated <- tibble::tibble(
    Kx = 2
  )

  before <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    P = 10,
    premium_unit = "annualized",
    benefit = 100,
    payment = 1,
    k = 1,
    type = "whole",
    annuity_type = "whole",
    premium_timing = "due",
    reserve_timing = "before_payment"
  )

  after <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    P = 10,
    premium_unit = "annualized",
    benefit = 100,
    payment = 1,
    k = 1,
    type = "whole",
    annuity_type = "whole",
    premium_timing = "due",
    reserve_timing = "after_payment"
  )

  expect_equal(before$Y_t, 3)
  expect_equal(after$Y_t, 2)
  expect_equal(before$future_pv_premiums, 30)
  expect_equal(after$future_pv_premiums, 20)
  expect_equal(before$L_t, 70)
  expect_equal(after$L_t, 80)
})


test_that("mc_reserve accepts standardized premium columns from mc_premium", {
  simulated <- tibble::tibble(
    Kx = c(1, 2),
    Tx = c(1.5, 2.5),
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12,
    premium_annualized = 120,
    premium_per_payment = 10,
    premium_annuity_scale = 120
  )

  result <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    benefit = 100,
    type = "whole",
    annuity_type = "whole"
  )

  expect_equal(unique(result$payments_per_year), 12L)
  expect_equal(unique(result$payment_per_payment), 1 / 12)
  expect_equal(unique(result$premium_annualized), 120)
  expect_equal(unique(result$premium_per_payment), 10)
  expect_equal(unique(result$premium_annuity_scale), 120)

  expect_equal(
    result$future_pv_premiums,
    result$premium_annuity_scale * result$Y_t
  )
})


test_that("mc_reserve legacy P remains an annuity-scale coefficient", {
  simulated <- tibble::tibble(
    Kx = 2,
    P = 10
  )

  result <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    benefit = 100,
    payment = 1,
    k = 1,
    type = "whole",
    annuity_type = "whole"
  )

  expect_equal(result$premium_annuity_scale, 10)
  expect_equal(result$premium_per_payment, 10)
  expect_equal(result$premium_annualized, 10)
  expect_equal(result$L_t, 70)
})


test_that("mc_reserve preserves missing lifetime information", {
  simulated <- tibble::tibble(
    Kx = c(2, NA_real_)
  )

  result <- mc_reserve(
    simulated,
    t = 0,
    i = 0.05,
    P = 10,
    premium_unit = "annualized",
    benefit = 100,
    payment = 1,
    k = 1,
    type = "whole",
    annuity_type = "whole",
    not_in_force = "zero"
  )

  expect_false(is.na(result$L_t[[1]]))
  expect_true(is.na(result$in_force[[2]]))
  expect_true(is.na(result$Z_t[[2]]))
  expect_true(is.na(result$Y_t[[2]]))
  expect_true(is.na(result$L_t[[2]]))
})


test_that("mc_reserve reports the premium payment timing correctly", {
  simulated <- tibble::tibble(
    Kx = 2
  )

  result <- mc_reserve(
    simulated,
    t = 0,
    i = 0,
    P = 10,
    premium_unit = "annualized",
    benefit = 100,
    type = "whole",
    annuity_type = "whole",
    premium_timing = "immediate"
  )

  expect_identical(result$payment_timing, "immediate")
  expect_identical(result$premium_timing, "immediate")
  expect_identical(result$timing, "end_of_year")
})


test_that("mc_reserve validates discrete terms and conversion frequency", {
  simulated <- tibble::tibble(
    Kx = c(1, 2),
    Tx = c(1.4, 2.6)
  )

  expect_error(
    mc_reserve(
      simulated,
      t = 0,
      i = 0.05,
      type = "term",
      n = 2.5,
      annuity_type = "temporary",
      timing = "end_of_year"
    ),
    "`n` must be a whole number"
  )

  expect_error(
    mc_reserve(
      simulated,
      t = 0,
      i = 0.06,
      i_type = "nominal_interest",
      m = 1.5,
      type = "whole",
      annuity_type = "whole"
    ),
    "`m`"
  )
})


test_that("mc_reserve validates premium schedule alignment", {
  simulated <- tibble::tibble(
    Kx = c(1, 2),
    Tx = c(1.4, 2.6)
  )

  expect_error(
    mc_reserve(
      simulated,
      t = 0,
      i = 0.05,
      type = "whole",
      annuity_type = "temporary",
      n = 2.3,
      k = 4,
      payment = 1 / 4
    ),
    "`n \\* k` must be an integer"
  )
})


test_that("mc_reserve detects inconsistent premium metadata", {
  simulated <- tibble::tibble(
    Kx = 2,
    Tx = 2.5,
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12,
    premium_annualized = 120,
    premium_per_payment = 11
  )

  expect_error(
    mc_reserve(
      simulated,
      t = 0,
      i = 0.05,
      benefit = 100,
      type = "whole",
      annuity_type = "whole"
    ),
    "inconsistent"
  )
})


test_that("mc_reserve deprecated aliases use exact conflict detection", {
  simulated <- tibble::tibble(
    Kx = 2
  )

  expect_error(
    mc_reserve(
      simulated,
      t = 0,
      duration = 0,
      i = 0.05
    ),
    "Provide only one of `t`"
  )

  expect_error(
    mc_reserve(
      simulated,
      i = 0.05,
      k = 1,
      payments_per_year = 1
    ),
    "Provide only one of `k`"
  )

  result <- mc_reserve(
    data = simulated,
    duration = 0,
    rate = 0.05,
    premium = 10,
    insurance = "whole_life",
    annuity = "whole_life"
  )

  expect_s3_class(result, "tbl_df")
})


test_that("mc_reserve requires an explicitly selected premium column to exist", {
  simulated <- tibble::tibble(
    Kx = 2
  )

  expect_error(
    mc_reserve(
      simulated,
      t = 0,
      i = 0.05,
      col_P = "missing_premium"
    ),
    "was not found"
  )
})
