test_that("premium_xy applies the equivalence principle for annual premiums", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  apv_benefits <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    n = 20,
    benefit = 100000,
    frac = "UDD"
  )

  apv_premium_annuity <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 10,
    k = 1,
    timing = "due",
    frac = "UDD"
  )

  expected <- apv_benefits / apv_premium_annuity

  observed <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("premium_xy distinguishes annualized and monthly premiums", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "last",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    frac = "UDD",
    output = "summary"
  )

  expect_equal(
    result$premium_annualized,
    12 * result$premium_per_payment,
    tolerance = 1e-12
  )

  expect_lt(abs(result$equivalence_residual), 1e-9)
})


test_that("premium_xy summary is compact", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "endowment",
    status = "joint",
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
})


test_that("piped and direct premium_xy calculations agree", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  direct <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
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
    lives = "joint",
    x = 60,
    y = 62,
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
    premium_xy(output = "summary")

  expect_equal(piped, direct, tolerance = 1e-12)
})


test_that("last-survivor contract infers status in the pipe", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  piped <- life_contract(
    lt = lt,
    lives = "last_survivor",
    x = 60,
    y = 62,
    i = 0.05
  ) |>
    add_insurance(
      type = "term",
      benefit = 100000,
      n = 20
    ) |>
    add_premium_schedule(
      k = 12,
      n_prem = 10
    ) |>
    premium_xy(output = "summary")

  direct <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "last",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    output = "summary"
  )

  expect_equal(piped, direct, tolerance = 1e-12)
})


test_that("premium_xy supports pure endowment benefits", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  pure <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "pure_endowment",
    status = "joint",
    benefit = 100000,
    n = 15,
    k = 1,
    n_prem = 10,
    output = "summary"
  )

  endowment_apv <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "endowment",
    status = "joint",
    benefit = 100000,
    n = 15,
    frac = "UDD"
  )

  term_apv <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 15,
    frac = "UDD"
  )

  expect_equal(
    pure$apv_benefits,
    endowment_apv - term_apv,
    tolerance = 1e-12
  )
})


test_that("premium_xy accepts fractional premium terms aligned with k", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 2.5,
    output = "summary"
  )

  expect_true(is.finite(result$premium_annualized))
  expect_true(is.finite(result$premium_per_payment))
})


test_that("premium_xy validates the actual premium-payment endpoint", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  valid <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 10,
    h = 5,
    k = 1,
    n_prem = 12,
    premium_start = "issue"
  )

  expect_true(is.finite(valid))

  expect_error(
    premium_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = 0.05,
      type = "term",
      status = "joint",
      benefit = 100000,
      n = 10,
      h = 5,
      k = 1,
      n_prem = 11,
      premium_start = "deferred"
    ),
    "must not extend beyond the end of coverage"
  )
})


test_that("premium_xy audit output is compact long format", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
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


test_that("deprecated tidy argument maps to premium_xy output levels", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  value_result <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    tidy = FALSE
  )

  summary_result <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    tidy = TRUE
  )

  expect_type(value_result, "double")
  expect_s3_class(summary_result, "tbl_df")
  expect_equal(ncol(summary_result), 6L)
})
