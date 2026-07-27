test_that("mc_premium returns annualized and per-payment premiums", {
  simulated <- tibble::tibble(
    pv_benefit = c(800, 1200),
    pv_annuity = c(8, 12),
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12
  )

  result <- mc_premium(simulated)

  expect_equal(result$premium_annuity_scale, c(100, 100))
  expect_equal(result$premium_annualized, c(100, 100))
  expect_equal(
    result$premium_per_payment,
    rep(100 / 12, 2),
    tolerance = 1e-12
  )
  expect_equal(result$P, result$premium_annuity_scale)
  expect_equal(result$premium, result$premium_annuity_scale)
  expect_lt(abs(result$mc_equivalence_residual[[1]]), 1e-12)
})


test_that("mc_premium is invariant to an equivalent annuity payment basis", {
  normalized <- tibble::tibble(
    pv_benefit = c(800, 1200),
    pv_annuity = c(8, 12),
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12
  )

  unit_per_payment <- tibble::tibble(
    pv_benefit = c(800, 1200),
    pv_annuity = c(96, 144),
    payment_per_payment = 1,
    payment_annualized = 12,
    payments_per_year = 12
  )

  result_normalized <- mc_premium(normalized)
  result_unit <- mc_premium(unit_per_payment)

  expect_equal(
    result_normalized$premium_annualized,
    result_unit$premium_annualized
  )

  expect_equal(
    result_normalized$premium_per_payment,
    result_unit$premium_per_payment
  )

  expect_equal(result_normalized$premium_annualized[[1]], 100)
  expect_equal(result_unit$premium_annuity_scale[[1]], 100 / 12)
})


test_that("mc_premium removes missing simulations as complete pairs", {
  simulated <- tibble::tibble(
    pv_benefit = c(1000, NA_real_, 600),
    pv_annuity = c(10, 20, 6)
  )

  result <- mc_premium(simulated)

  expect_equal(result$premium_annualized, rep(100, 3))

  expect_error(
    mc_premium(simulated, na_rm = FALSE),
    "Missing simulated values"
  )
})


test_that("mc_premium computes premiums by group", {
  simulated <- tibble::tibble(
    age = c(40, 40, 50, 50),
    pv_benefit = c(800, 1200, 900, 1500),
    pv_annuity = c(8, 12, 6, 10),
    payment_per_payment = 1 / 4,
    payment_annualized = 1,
    payments_per_year = 4
  )

  result <- mc_premium(simulated, by = "age")

  premium_by_age <- result |>
    dplyr::distinct(age, premium_annualized) |>
    dplyr::arrange(age)

  expect_equal(
    premium_by_age$premium_annualized,
    c(100, 150)
  )
})


test_that("mc_premium works with an existing dplyr grouping", {
  simulated <- tibble::tibble(
    age = c(40, 40, 50, 50),
    pv_benefit = c(800, 1200, 900, 1500),
    pv_annuity = c(8, 12, 6, 10)
  )

  result <- simulated |>
    dplyr::group_by(age) |>
    mc_premium()

  expect_false(dplyr::is_grouped_df(result))

  premium_by_age <- result |>
    dplyr::distinct(age, premium_annualized) |>
    dplyr::arrange(age)

  expect_equal(
    premium_by_age$premium_annualized,
    c(100, 150)
  )
})


test_that("mc_premium and mc_loss satisfy Monte Carlo equivalence", {
  simulated <- tibble::tibble(
    pv_benefit = c(800, 1200, 1000),
    pv_annuity = c(8, 12, 10),
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12
  )

  result <- simulated |>
    mc_premium() |>
    mc_loss()

  expect_lt(abs(mean(result$L)), 1e-12)
  expect_equal(
    result$premium_annualized,
    rep(100, 3)
  )
})


test_that("mc_premium supports legacy annual inputs", {
  simulated <- tibble::tibble(
    Z = c(800, 1200),
    Y = c(8, 12)
  )

  result <- mc_premium(
    data = simulated,
    benefit_col = "Z",
    annuity_col = "Y",
    premium_col = "premium"
  )

  expect_equal(result$premium, c(100, 100))
  expect_equal(result$premium_annuity_scale, c(100, 100))
  expect_equal(result$premium_annualized, c(100, 100))
  expect_equal(result$premium_per_payment, c(100, 100))
})


test_that("mc_premium deprecated aliases use exact conflict detection", {
  simulated <- tibble::tibble(
    Z = c(800, 1200),
    Y = c(8, 12)
  )

  expect_error(
    mc_premium(
      simulated,
      col_Z = "Z",
      benefit_col = "Z",
      col_Y = "Y"
    ),
    "Provide only one of `col_Z`"
  )
})


test_that("mc_premium requires fractional annuity payment metadata", {
  simulated <- tibble::tibble(
    pv_benefit = c(800, 1200),
    pv_annuity = c(8, 12),
    payments_per_year = 12
  )

  expect_error(
    mc_premium(simulated),
    "cannot be inferred"
  )

  result <- mc_premium(
    simulated,
    annuity_payment = 1 / 12
  )

  expect_equal(result$premium_annualized, c(100, 100))
})


test_that("mc_premium detects inconsistent payment metadata", {
  simulated <- tibble::tibble(
    pv_benefit = c(800, 1200),
    pv_annuity = c(8, 12),
    payment_per_payment = 1 / 12,
    payment_annualized = 2,
    payments_per_year = 12
  )

  expect_error(
    mc_premium(simulated),
    "payment_annualized"
  )
})


test_that("mc_premium rejects standardized names for col_P", {
  simulated <- tibble::tibble(
    pv_benefit = c(800, 1200),
    pv_annuity = c(8, 12)
  )

  expect_error(
    mc_premium(
      simulated,
      col_P = "premium_annualized"
    ),
    "standardized premium column"
  )
})
