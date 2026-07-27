test_that("contract components support joint-life pipelines", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )

  contract <- life_contract(
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
      premium_start = "issue"
    )

  expect_identical(contract$insurance$type, "term")
  expect_identical(contract$insurance$status, "joint")
  expect_identical(contract$insurance$benefit, 100000)
  expect_identical(contract$premium_schedule$k, 12L)
  expect_identical(contract$premium_schedule$n_prem, 10)
})


test_that("contract components support last-survivor pipelines", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )

  contract <- life_contract(
    lt = lt,
    lives = "last_survivor",
    x = 60,
    y = 62,
    i = 0.05
  ) |>
    add_insurance(
      type = "pure_endowment",
      benefit = 50000,
      n = 15
    ) |>
    add_premium_schedule(
      k = 4,
      n_prem = 10
    )

  expect_identical(contract$insurance$type, "pure_endowment")
  expect_identical(contract$insurance$status, "last")
  expect_identical(contract$premium_schedule$k, 4L)
})


test_that("contract components reject incompatible insurance types", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )

  two_life <- life_contract(
    lt = lt,
    lives = "joint",
    x = 60,
    y = 62,
    i = 0.05
  )

  single <- life_contract(
    lt = lt,
    lives = "single",
    x = 60,
    i = 0.05
  )

  expect_error(
    add_insurance(
      two_life,
      type = "variable_k",
      benefit = rep(1000, 12),
      n = 1,
      k = 12
    ),
    "not supported for a two-life contract"
  )

  expect_error(
    add_insurance(
      single,
      type = "pure_endowment",
      benefit = 1000,
      n = 10
    ),
    "not supported for a single-life contract"
  )
})


test_that("two-life premium schedules accept fractional terms aligned with k", {
  lt <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.012 * (0:60)^1.35))
  )

  contract <- life_contract(
    lt = lt,
    lives = "joint",
    x = 60,
    y = 62,
    i = 0.05
  ) |>
    add_premium_schedule(
      k = 12,
      n_prem = 2.5
    )

  expect_identical(contract$premium_schedule$n_prem, 2.5)
})
