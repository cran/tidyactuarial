test_that("reserve_x prospective reserve follows the future-loss formula", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  P <- premium_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10
  )

  expected_benefits <- insurance_x(
    lt = lt,
    x = 45,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 15
  )

  expected_premiums <- P * annuity_x(
    lt = lt,
    x = 45,
    i = 0.05,
    n = 5,
    k = 1,
    timing = "due"
  )

  observed <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = 5,
    output = "value"
  )

  expect_equal(
    unname(observed),
    expected_benefits - expected_premiums,
    tolerance = 1e-10
  )
})


test_that("reserve_x uses the annualized premium for true monthly reserves", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  premium <- premium_x(
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

  expected_benefits <- insurance_x(
    lt = lt,
    x = 45,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 15
  )

  expected_factor <- annuity_x(
    lt = lt,
    x = 45,
    i = 0.05,
    n = 5,
    k = 12,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  expected <- expected_benefits -
    premium$premium_annualized * expected_factor

  observed <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    frac = "UDD",
    t = 5,
    method = "prospective",
    output = "value"
  )

  wrong_using_installment <- expected_benefits -
    premium$premium_per_payment * expected_factor

  expect_equal(unname(observed), expected, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(observed), wrong_using_installment)))
})


test_that("reserve_x net reserve at issue is zero", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "endowment",
    benefit = 50000,
    n = 20,
    k = 12,
    n_prem = 10,
    t = 0,
    output = "value"
  )

  expect_lt(abs(unname(result)), 1e-8)
})


test_that("piped and direct reserve_x calculations agree", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  direct <- reserve_x(
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
    t = c(0, 5, 10, 15, 20),
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
    reserve_x(
      t = c(0, 5, 10, 15, 20),
      output = "summary"
    )

  expect_equal(piped, direct, tolerance = 1e-10)
})


test_that("annual recursive and prospective reserves agree", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  durations <- 0:20

  prospective <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "endowment",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = durations,
    method = "prospective",
    output = "value"
  )

  recursive <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "endowment",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = durations,
    method = "recursive",
    output = "value"
  )

  expect_equal(recursive, prospective, tolerance = 1e-7)
})


test_that("endowment reserve at maturity equals the benefit", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "endowment",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = 20,
    output = "value"
  )

  expect_equal(unname(result), 100000, tolerance = 1e-8)
})


test_that("reserve_x summary remains compact and identifies premium units", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    t = c(0, 5, 10),
    output = "summary"
  )

  expect_identical(
    names(result),
    c(
      "t",
      "age",
      "reserve",
      "premium_annualized",
      "premium_per_payment",
      "method"
    )
  )

  expect_equal(ncol(result), 6L)

  expect_equal(
    result$premium_annualized,
    12 * result$premium_per_payment,
    tolerance = 1e-12
  )
})


test_that("reserve_x audit output exposes prospective components", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_x(
    lt = lt,
    x = 40,
    i = 0.05,
    type = "term",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    t = c(0, 5),
    output = "audit"
  )

  expect_identical(
    names(result),
    c("t", "age", "component", "value", "unit")
  )

  expect_true(
    all(
      c(
        "apv_future_benefits",
        "apv_future_premiums",
        "reserve",
        "premium_annualized",
        "premium_per_payment"
      ) %in% result$component
    )
  )
})


test_that("recursive reserves reject k-thly premium schedules", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  expect_error(
    reserve_x(
      lt = lt,
      x = 40,
      i = 0.05,
      type = "term",
      benefit = 100000,
      n = 20,
      k = 12,
      n_prem = 10,
      method = "recursive"
    ),
    "recursive method currently supports annual premiums only"
  )
})
