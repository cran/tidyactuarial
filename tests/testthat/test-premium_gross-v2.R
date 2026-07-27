test_that("premium_gross preserves the previous annual formula", {
  prem <- tibble::tibble(
    premium_annualized = 1200,
    premium_per_payment = 1200,
    payments_per_year = 1,
    apv_premium_annuity = 12.5
  )

  expected <- (1200 + 50) / ((1 - 0.05) - 0.5 / 12.5)

  observed <- premium_gross(
    prem,
    alpha = 0.5,
    beta = 0.05,
    gamma = 50
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("premium_gross correctly scales k-thly expense units", {
  prem <- tibble::tibble(
    premium_annualized = 1200,
    premium_per_payment = 100,
    payments_per_year = 12,
    apv_premium_annuity = 10
  )

  expected <- (1200 + 12 * 20) /
    ((1 - 0.05) - 0.5 / (12 * 10))

  result <- premium_gross(
    prem,
    alpha = 0.5,
    beta = 0.05,
    gamma = 20,
    output = "summary"
  )

  expect_equal(
    result$gross_premium_annualized,
    expected,
    tolerance = 1e-12
  )

  expect_equal(
    result$gross_premium_per_payment,
    expected / 12,
    tolerance = 1e-12
  )

  expect_lt(abs(result$equivalence_residual), 1e-9)
})


test_that("premium_gross expense APVs satisfy extended equivalence", {
  prem <- tibble::tibble(
    premium_annualized = 2400,
    premium_per_payment = 200,
    payments_per_year = 12,
    apv_premium_annuity = 8.75
  )

  audit <- premium_gross(
    prem,
    alpha = 1.25,
    beta = 0.08,
    gamma = 15,
    output = "audit"
  )

  value_of <- function(name) {
    audit$value[audit$component == name]
  }

  expect_equal(
    value_of("apv_gross_premiums"),
    value_of("apv_benefits") +
      value_of("apv_total_expenses"),
    tolerance = 1e-9
  )

  expect_equal(
    value_of("apv_initial_expense"),
    1.25 * value_of("gross_premium_per_payment"),
    tolerance = 1e-12
  )

  expect_equal(
    value_of("apv_fixed_per_payment_expense"),
    15 * 12 * 8.75,
    tolerance = 1e-12
  )
})


test_that("premium_gross accepts a premium_x summary", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  net <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    output = "summary"
  )

  gross <- premium_gross(
    net,
    alpha = 0.5,
    beta = 0.05,
    gamma = 10,
    output = "summary"
  )

  expect_s3_class(gross, "tbl_df")
  expect_equal(ncol(gross), 6L)
  expect_gt(
    gross$gross_premium_annualized,
    gross$net_premium_annualized
  )
  expect_equal(
    gross$gross_premium_annualized,
    12 * gross$gross_premium_per_payment,
    tolerance = 1e-12
  )
})


test_that("premium_gross accepts a premium_xy summary", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  net <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 4,
    n_prem = 10,
    output = "summary"
  )

  gross <- premium_gross(
    net,
    alpha = 1,
    beta = 0.04,
    gamma = 12,
    output = "summary"
  )

  expect_s3_class(gross, "tbl_df")
  expect_equal(gross$payments_per_year, 4L)
  expect_lt(abs(gross$equivalence_residual), 1e-8)
})


test_that("premium_gross supports a pipe from life_contract", {
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
      n_prem = 10
    ) |>
    premium_x(output = "summary") |>
    premium_gross(
      alpha = 0.5,
      beta = 0.05,
      gamma = 10,
      output = "summary"
    )

  expect_s3_class(result, "tbl_df")
  expect_equal(result$payments_per_year, 12L)
})


test_that("premium_gross summary is compact", {
  prem <- tibble::tibble(
    premium_annualized = 1200,
    premium_per_payment = 100,
    payments_per_year = 12,
    apv_premium_annuity = 10
  )

  result <- premium_gross(
    prem,
    alpha = 0.5,
    beta = 0.05,
    gamma = 20,
    output = "summary"
  )

  expect_identical(
    names(result),
    c(
      "gross_premium_annualized",
      "gross_premium_per_payment",
      "net_premium_annualized",
      "loading_annualized",
      "payments_per_year",
      "equivalence_residual"
    )
  )

  expect_equal(ncol(result), 6L)
})


test_that("premium_gross accepts legacy annual tables", {
  prem <- tibble::tibble(
    P = 1200,
    a_premiums = 12.5
  )

  observed <- premium_gross(
    prem,
    alpha = 0.5,
    beta = 0.05,
    gamma = 50
  )

  expected <- (1200 + 50) / ((1 - 0.05) - 0.5 / 12.5)

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("premium_gross rejects ambiguous legacy subannual premiums", {
  prem <- tibble::tibble(
    P = 100,
    payments_per_year = 12,
    a_premiums = 10
  )

  expect_error(
    premium_gross(prem),
    "Legacy premium columns are ambiguous"
  )
})


test_that("premium_gross detects inconsistent premium units", {
  prem <- tibble::tibble(
    premium_annualized = 1200,
    premium_per_payment = 90,
    payments_per_year = 12,
    apv_premium_annuity = 10
  )

  expect_error(
    premium_gross(prem),
    "inconsistent"
  )
})


test_that("premium_gross rejects a nonpositive denominator", {
  prem <- tibble::tibble(
    premium_annualized = 1200,
    premium_per_payment = 1200,
    payments_per_year = 1,
    apv_premium_annuity = 1
  )

  expect_error(
    premium_gross(
      prem,
      alpha = 1,
      beta = 0,
      gamma = 0
    ),
    "denominator nonpositive"
  )
})


test_that("premium_gross allows a zero net premium with fixed expenses", {
  prem <- tibble::tibble(
    premium_annualized = 0,
    premium_per_payment = 0,
    payments_per_year = 12,
    apv_premium_annuity = 10
  )

  result <- premium_gross(
    prem,
    gamma = 5,
    output = "summary"
  )

  expect_gt(result$gross_premium_annualized, 0)
  expect_equal(
    result$gross_premium_annualized,
    12 * 5,
    tolerance = 1e-12
  )
})


test_that("deprecated tidy maps to premium_gross output levels", {
  prem <- tibble::tibble(
    premium_annualized = 1200,
    premium_per_payment = 1200,
    payments_per_year = 1,
    apv_premium_annuity = 12.5
  )

  value_result <- premium_gross(prem, tidy = FALSE)
  summary_result <- premium_gross(prem, tidy = TRUE)

  expect_type(value_result, "double")
  expect_s3_class(summary_result, "tbl_df")
  expect_equal(ncol(summary_result), 6L)
})
