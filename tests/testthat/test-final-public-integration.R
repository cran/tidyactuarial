make_lt_final_integration <- function() {
  lifetable(
    x = 60:65,
    lx = c(100000, 90000, 75000, 55000, 30000, 0),
    close = TRUE,
    frac = "UDD"
  )
}

test_that("A_xj remains a numerical alias of insurance_xj", {
  md <- md_table(
    tibble::tibble(
      x = 60:63,
      q_death = c(0.10, 0.20, 0.30, 1.00),
      q_other = c(0.10, 0.10, 0.10, 0.00)
    ),
    radix = 1000,
    close = TRUE
  )

  current <- insurance_xj(
    md = md,
    x = 60,
    i = 0.05,
    cause = "q_death",
    type = "term",
    n = 3
  )

  alias <- A_xj(
    md = md,
    x = 60,
    i = 0.05,
    cause = "q_death",
    type = "term",
    n = 3
  )

  expect_equal(alias, current, tolerance = 1e-12)
})

test_that("simulate_annuity_x summary output equals summary_mc of the simulations", {
  lt <- make_lt_final_integration()

  sim <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 4,
    k = 2,
    timing = "due",
    payment = 0.5,
    n_sim = 100,
    seed = 31415,
    output = "simulations"
  )

  expected <- summary_mc(
    sim,
    col_value = "present_value"
  )

  actual <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 4,
    k = 2,
    timing = "due",
    payment = 0.5,
    n_sim = 100,
    seed = 31415,
    output = "summary"
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("simulate_insurance_x summary output equals summary_mc of the simulations", {
  lt <- make_lt_final_integration()

  sim <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0.05,
    type = "term",
    n = 4,
    benefit = 100000,
    n_sim = 100,
    seed = 27182,
    output = "simulations"
  )

  expected <- summary_mc(
    sim,
    col_value = "present_value"
  )

  actual <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0.05,
    type = "term",
    n = 4,
    benefit = 100000,
    n_sim = 100,
    seed = 27182,
    output = "summary"
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("fixed seeds make the single-life simulation interfaces reproducible", {
  lt <- make_lt_final_integration()

  a1 <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0.04,
    n = 3,
    n_sim = 50,
    seed = 123
  )

  a2 <- simulate_annuity_x(
    lt = lt,
    x = 60,
    i = 0.04,
    n = 3,
    n_sim = 50,
    seed = 123
  )

  expect_equal(a1, a2)

  i1 <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0.04,
    type = "endowment",
    n = 3,
    n_sim = 50,
    seed = 456
  )

  i2 <- simulate_insurance_x(
    lt = lt,
    x = 60,
    i = 0.04,
    type = "endowment",
    n = 3,
    n_sim = 50,
    seed = 456
  )

  expect_equal(i1, i2)
})
