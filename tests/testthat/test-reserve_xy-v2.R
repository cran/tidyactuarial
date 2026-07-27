test_that("reserve_xy follows the prospective future-loss formula", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  P <- premium_xy(
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

  expected <- insurance_xy(
    lt = lt,
    x = 65,
    y = 67,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 15,
    frac = "UDD"
  ) -
    P * annuity_xy(
      lt = lt,
      x = 65,
      y = 67,
      i = 0.05,
      status = "joint",
      n = 5,
      k = 1,
      timing = "due",
      frac = "UDD"
    )

  observed <- reserve_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = 5,
    output = "value"
  )

  expect_equal(unname(observed), expected, tolerance = 1e-10)
})


test_that("reserve_xy uses the annualized premium for monthly reserves", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  premium <- premium_xy(
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

  benefit_apv <- insurance_xy(
    lt = lt,
    x = 65,
    y = 67,
    i = 0.05,
    type = "term",
    status = "last",
    benefit = 100000,
    n = 15,
    frac = "UDD"
  )

  premium_factor <- annuity_xy(
    lt = lt,
    x = 65,
    y = 67,
    i = 0.05,
    status = "last",
    n = 5,
    k = 12,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  expected <- benefit_apv -
    premium$premium_annualized * premium_factor

  wrong <- benefit_apv -
    premium$premium_per_payment * premium_factor

  observed <- reserve_xy(
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
    t = 5,
    output = "value"
  )

  expect_equal(unname(observed), expected, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(observed), wrong)))
})


test_that("reserve_xy net reserve at issue is zero", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "endowment",
    status = "joint",
    benefit = 50000,
    n = 20,
    k = 12,
    n_prem = 10,
    t = 0,
    output = "value"
  )

  expect_lt(abs(unname(result)), 1e-8)
})


test_that("direct and piped reserve_xy calculations agree", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  direct <- reserve_xy(
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
    frac = "UDD",
    t = c(0, 5, 10, 15, 20),
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
      frac = "UDD"
    ) |>
    add_premium_schedule(
      k = 12,
      n_prem = 10,
      frac = "UDD"
    ) |>
    reserve_xy(
      t = c(0, 5, 10, 15, 20),
      output = "summary"
    )

  expect_equal(piped, direct, tolerance = 1e-10)
})


test_that("last-survivor status is inferred from the contract", {
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
    reserve_xy(
      t = c(0, 5, 10),
      output = "value"
    )

  direct <- reserve_xy(
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
    t = c(0, 5, 10),
    output = "value"
  )

  expect_equal(piped, direct, tolerance = 1e-10)
})


test_that("annual recursive and prospective reserves agree", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  durations <- 0:20

  prospective <- reserve_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "endowment",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = durations,
    method = "prospective",
    output = "value"
  )

  recursive <- reserve_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "endowment",
    status = "joint",
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


test_that("a supplied non-net premium is preserved at issue", {
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
    k = 1,
    n_prem = 10
  )

  result <- reserve_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    P = 1.10 * net,
    k = 1,
    n_prem = 10,
    t = 0,
    output = "value"
  )

  expect_gt(abs(unname(result)), 1e-8)
})


test_that("pure endowment reserve at maturity equals the benefit", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_xy(
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
    t = 15,
    output = "value"
  )

  expect_equal(unname(result), 100000, tolerance = 1e-8)
})


test_that("reserve_xy summary is compact and labels premium units", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_xy(
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
    t = c(0, 5, 10),
    output = "summary"
  )

  expect_identical(
    names(result),
    c(
      "t", "age_x", "age_y", "reserve",
      "premium_annualized", "premium_per_payment"
    )
  )
  expect_equal(ncol(result), 6L)
  expect_equal(
    result$premium_annualized,
    12 * result$premium_per_payment,
    tolerance = 1e-12
  )
})


test_that("reserve_xy audit exposes the prospective components", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  result <- reserve_xy(
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
    t = c(0, 5),
    output = "audit"
  )

  expect_identical(
    names(result),
    c("t", "age_x", "age_y", "component", "value", "unit")
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


test_that("recursive reserve_xy rejects k-thly and deferred structures", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  expect_error(
    reserve_xy(
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
      method = "recursive"
    ),
    "recursive method currently requires"
  )

  expect_error(
    reserve_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = 0.05,
      type = "term",
      status = "joint",
      benefit = 100000,
      n = 20,
      h = 2,
      k = 1,
      n_prem = 10,
      method = "recursive"
    ),
    "recursive method currently requires"
  )
})
