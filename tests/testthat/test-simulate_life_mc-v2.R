test_that("simulate_annuity_x uses per-payment amounts explicitly", {
  lt <- data.frame(
    x = 60:61,
    lx = c(100, 0)
  )

  result <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0,
    n = 1,
    k = 12,
    payment = 1 / 12,
    frac = "udd",
    timing = "due",
    n_sim = 200,
    seed = 123
  )

  expect_equal(result$payment_per_payment, rep(1 / 12, 200))
  expect_equal(result$payment_annualized, rep(1, 200))
  expect_equal(
    result$present_value,
    result$n_payments / 12,
    tolerance = 1e-12
  )
  expect_equal(result$pv_annuity, result$present_value)
})


test_that("simulate_annuity_x uses complete lifetimes for fractional payments", {
  lt <- data.frame(
    x = 60:61,
    lx = c(100, 0)
  )

  result <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0,
    n = 1,
    k = 12,
    payment = 1 / 12,
    frac = "udd",
    timing = "due",
    n_sim = 500,
    seed = 321
  )

  expect_true(all(result$Kx == 0L))
  expect_true(all(result$Tx >= 0 & result$Tx < 1))
  expect_gt(length(unique(result$n_payments)), 1L)
  expect_true(all(result$n_payments >= 1L))
  expect_true(all(result$n_payments <= 12L))
})


test_that("simulate_annuity_x preserves annual due and immediate conventions", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 100, 0)
  )

  due <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0,
    k = 1,
    timing = "due",
    payment = 10,
    n_sim = 20,
    seed = 1
  )

  immediate <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0,
    k = 1,
    timing = "immediate",
    payment = 10,
    n_sim = 20,
    seed = 1
  )

  expect_equal(due$Kx, rep(1L, 20))
  expect_equal(due$n_payments, rep(2L, 20))
  expect_equal(due$present_value, rep(20, 20))

  expect_equal(immediate$n_payments, rep(1L, 20))
  expect_equal(immediate$present_value, rep(10, 20))
})


test_that("fractional lifetime assumptions produce valid complete lifetimes", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 60, 0)
  )

  for (assumption in c("udd", "cml", "balducci")) {
    result <- simulate_annuity_x(
      lt = lt,
      x = 60,
      i = 0.05,
      k = 4,
      payment = 1 / 4,
      frac = assumption,
      n_sim = 200,
      seed = 100
    )

    expect_true(
      all(
        result$Tx >= result$Kx &
          result$Tx < result$Kx + 1
      )
    )
  }
})


test_that("simulate_insurance_x produces coherent benefit outputs", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 50, 0)
  )

  whole <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0,
    type = "whole",
    benefit = 1000,
    n_sim = 200,
    seed = 123
  )

  endowment <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0,
    type = "endowment",
    n = 1,
    benefit = 1000,
    n_sim = 200,
    seed = 123
  )

  expect_equal(whole$benefit_indicator, rep(1, 200))
  expect_equal(whole$present_value, rep(1000, 200))
  expect_equal(whole$pv_benefit, whole$present_value)

  expect_equal(endowment$benefit_indicator, rep(1, 200))
  expect_equal(endowment$present_value, rep(1000, 200))
  expect_equal(endowment$pv_benefit, endowment$present_value)
})


test_that("simulate_insurance_x term coverage follows the simulated death year", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 50, 0)
  )

  result <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0,
    type = "term",
    n = 1,
    benefit = 1000,
    n_sim = 5000,
    seed = 42
  )

  expect_equal(
    result$benefit_indicator,
    as.numeric(result$Kx < 1 & result$died_within_horizon)
  )

  expect_equal(
    result$present_value,
    1000 * result$benefit_indicator
  )

  expect_equal(
    mean(result$benefit_indicator),
    0.5,
    tolerance = 0.03
  )
})


test_that("life-table simulation rejects increasing lx", {
  invalid_lt <- data.frame(
    x = 60:62,
    lx = c(100, 110, 0)
  )

  expect_error(
    simulate_annuity_x(
      lt = invalid_lt,
      x = 60,
      i = 0.05,
      n_sim = 10
    ),
    "nonincreasing"
  )
})


test_that("seeded simulations restore the caller RNG state", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 50, 0)
  )

  set.seed(999)
  state_before <- .Random.seed

  invisible(
    simulate_annuity_x(
      lt = lt,
      x = 60,
      i = 0.05,
      n_sim = 20,
      seed = 123
    )
  )

  expect_identical(.Random.seed, state_before)
})


test_that("deprecated arguments use exact conflict detection", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 50, 0)
  )

  expect_error(
    simulate_annuity_x(
      lt = lt,
      mortality_table = lt,
      x = 60,
      i = 0.05,
      n_sim = 10
    ),
    "Provide only one of `lt`"
  )

  expect_error(
    simulate_annuity_x(
      lt = lt,
      x = 60,
      i = 0.05,
      k = 1,
      payments_per_year = 1,
      n_sim = 10
    ),
    "Provide only one of `k`"
  )

  expect_error(
    simulate_insurance_x(
      lt = lt,
      x = 60,
      i = 0.05,
      type = "whole",
      insurance_type = "whole",
      n_sim = 10
    ),
    "Provide only one of `type`"
  )
})


test_that("deprecated argument names remain usable", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 50, 0)
  )

  annuity <- simulate_annuity_x(
    mortality_table = lt,
    age = 60,
    rate = 0.05,
    rate_type = "effective",
    term_years = 1,
    payments_per_year = 1,
    n_sim = 10,
    seed = 1
  )

  insurance <- simulate_insurance_x(
    mortality_table = lt,
    age = 60,
    rate = 0.05,
    rate_type = "effective",
    insurance_type = "term",
    term_years = 1,
    n_sim = 10,
    seed = 1
  )

  expect_s3_class(annuity, "tbl_df")
  expect_s3_class(insurance, "tbl_df")
})


test_that(".annuity_factor_count does not hide a one-over-k normalization", {
  factor <- tidyactuarial:::.annuity_factor_count(
    n_payments = 12,
    i_effective = 0,
    k = 12,
    timing = "due"
  )

  expect_equal(factor, 12)
})


test_that("simulate_life_mc validates seeds and unnamed dots", {
  lt <- data.frame(
    x = 60:62,
    lx = c(100, 50, 0)
  )

  expect_error(
    simulate_annuity_x(
      lt = lt,
      x = 60,
      i = 0.05,
      seed = 1.5
    ),
    "nonnegative integer"
  )

  expect_error(
    tidyactuarial:::.life_mc_collect_old_args(
      dots = list(1),
      allowed = character()
    ),
    "must be named"
  )
})
