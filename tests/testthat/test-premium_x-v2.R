test_that("premium_x applies the equation of equivalence for annual premiums", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  apv_benefits <- insurance_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20
  )

  apv_premium_annuity <- annuity_x(
    lt = lt,
    x = 40,
    i = 0.05,
    n = 10,
    k = 1,
    timing = "due"
  )

  expected <- apv_benefits / apv_premium_annuity

  observed <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("premium_x distinguishes annualized and monthly premiums", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    frac = "UDD",
    output = "summary"
  )

  expected_annuity <- annuity_x(
    lt = lt,
    x = 40,
    i = 0.05,
    n = 10,
    k = 12,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  expected_benefits <- insurance_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20
  )

  expected_annualized <- expected_benefits / expected_annuity

  expect_equal(
    result$premium_annualized,
    expected_annualized,
    tolerance = 1e-12
  )

  expect_equal(
    result$premium_per_payment,
    expected_annualized / 12,
    tolerance = 1e-12
  )

  expect_equal(
    result$premium_annualized,
    12 * result$premium_per_payment,
    tolerance = 1e-12
  )
})


test_that("premium_x summary is compact and closes the equivalence equation", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "endowment",
    benefit = 50000,
    n = 15,
    k = 4,
    n_prem = 10,
    output = "summary"
  )

  expect_identical(
    names(result),
    c(
      "premium_annualized",
      "premium_per_payment",
      "payments_per_year",
      "apv_benefits",
      "apv_premium_annuity",
      "equivalence_residual"
    )
  )

  expect_equal(ncol(result), 6L)
  expect_lt(abs(result$equivalence_residual), 1e-9)
})


test_that("piped and direct premium_x calculations agree", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  direct <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    h = 0,
    k = 12,
    n_prem = 10,
    timing = "due",
    premium_start = "issue",
    frac = "UDD",
    output = "summary"
  )

  piped <- life_contract(
    lt = lt,
    lives = "single",
    x = 40,
    i = 0.05
  ) |>
    add_insurance(
      type = "term",
      benefit = 100000,
      n = 20,
      h = 0,
      frac = "UDD"
    ) |>
    add_premium_schedule(
      k = 12,
      n_prem = 10,
      timing = "due",
      premium_start = "issue",
      frac = "UDD"
    ) |>
    premium_x(output = "summary")

  expect_equal(piped, direct, tolerance = 1e-12)
})


test_that("premium_x supports fractional premium-paying terms", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- life_contract(
    lt = lt,
    lives = "single",
    x = 40,
    i = 0.05
  ) |>
    add_insurance(
      type = "term",
      benefit = 100000,
      n = 20
    ) |>
    add_premium_schedule(
      k = 12,
      n_prem = 2.5,
      timing = "due"
    ) |>
    premium_x(output = "summary")

  expect_true(is.finite(result$premium_annualized))
  expect_true(is.finite(result$premium_per_payment))
  expect_equal(
    result$premium_annualized,
    12 * result$premium_per_payment,
    tolerance = 1e-12
  )
})


test_that("premium_x audit output is compact long format", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    output = "audit"
  )

  expect_identical(names(result), c("component", "value", "unit"))
  expect_true(
    all(
      c(
        "premium_annualized",
        "premium_per_payment",
        "equivalence_residual"
      ) %in% result$component
    )
  )
})


test_that("deprecated tidy argument maps to the new output levels", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  value_result <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    tidy = FALSE
  )

  summary_result <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    tidy = TRUE
  )

  expect_type(value_result, "double")
  expect_s3_class(summary_result, "tbl_df")
  expect_equal(ncol(summary_result), 6L)
})
