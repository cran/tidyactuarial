make_md_insurance_test <- function() {
  md_table(
    tibble::tibble(
      x = 30:33,
      q_death = c(0.10, 0.20, 0.30, 1.00),
      q_disability = c(0.20, 0.10, 0.10, 0.00)
    ),
    radix = 1000,
    close = TRUE
  )
}

test_that("insurance_xj matches the direct cause-specific term-insurance sum", {
  md <- make_md_insurance_test()
  i <- 0.05
  v <- 1 / (1 + i)

  # Death in years 1, 2, 3:
  # q30^d + p30^(tau) q31^d + p30^(tau)p31^(tau) q32^d
  expected <- (
    v^1 * 0.10 +
    v^2 * 0.70 * 0.20 +
    v^3 * 0.70 * 0.70 * 0.30
  )

  expect_equal(
    insurance_xj(
      md = md,
      x = 30,
      i = i,
      cause = "q_death",
      type = "term",
      n = 3
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("insurance_xj applies deferment at the correct policy years", {
  md <- make_md_insurance_test()
  i <- 0.05
  v <- 1 / (1 + i)

  # Two-year term deferred one year: benefits can arise in years 2 and 3.
  expected <- (
    v^2 * 0.70 * 0.20 +
    v^3 * 0.70 * 0.70 * 0.30
  )

  expect_equal(
    insurance_xj(
      md = md,
      x = 30,
      i = i,
      cause = "q_death",
      type = "term",
      n = 2,
      h = 1
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("insurance_xj whole life spans the available multiple-decrement table", {
  md <- make_md_insurance_test()
  i <- 0.04
  v <- 1 / (1 + i)

  expected <- (
    v^1 * 0.10 +
    v^2 * 0.70 * 0.20 +
    v^3 * 0.70 * 0.70 * 0.30 +
    v^4 * 0.70 * 0.70 * 0.60 * 1.00
  )

  expect_equal(
    insurance_xj(
      md = md,
      x = 30,
      i = i,
      cause = "q_death",
      type = "whole"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("cause-specific insurance values add across mutually exclusive causes", {
  md <- make_md_insurance_test()
  i <- 0.03
  v <- 1 / (1 + i)

  death <- insurance_xj(
    md = md,
    x = 30,
    i = i,
    cause = "q_death",
    type = "term",
    n = 3
  )

  disability <- insurance_xj(
    md = md,
    x = 30,
    i = i,
    cause = "q_disability",
    type = "term",
    n = 3
  )

  expected_any_decrement <- (
    v^1 * 0.30 +
    v^2 * 0.70 * 0.30 +
    v^3 * 0.70 * 0.70 * 0.40
  )

  expect_equal(
    death + disability,
    expected_any_decrement,
    tolerance = 1e-12
  )
})

test_that("insurance_xj scales linearly with benefit", {
  md <- make_md_insurance_test()

  unit <- insurance_xj(
    md = md,
    x = 30,
    i = 0.05,
    cause = "q_death",
    type = "term",
    n = 3
  )

  scaled <- insurance_xj(
    md = md,
    x = 30,
    i = 0.05,
    cause = "q_death",
    type = "term",
    n = 3,
    benefit = 25000
  )

  expect_equal(scaled, 25000 * unit, tolerance = 1e-10)
})

test_that("insurance_xj supports equivalent nominal interest inputs", {
  md <- make_md_insurance_test()

  i_nominal <- 0.06
  m <- 12
  i_eff <- (1 + i_nominal / m)^m - 1

  from_nominal <- insurance_xj(
    md = md,
    x = 30,
    i = i_nominal,
    i_type = "nominal_interest",
    m = m,
    cause = "q_death",
    type = "term",
    n = 3
  )

  from_effective <- insurance_xj(
    md = md,
    x = 30,
    i = i_eff,
    i_type = "effective",
    cause = "q_death",
    type = "term",
    n = 3
  )

  expect_equal(from_nominal, from_effective, tolerance = 1e-12)
})

test_that("insurance_xj vectorizes and returns its audit tibble", {
  md <- make_md_insurance_test()

  out <- insurance_xj(
    md = md,
    x = c(30, 31),
    i = c(0.04, 0.05),
    cause = "q_death",
    type = "term",
    n = c(2, 2),
    benefit = c(1, 2),
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 2L)
  expect_named(
    out,
    c(
      "x", "i", "i_type", "m", "h", "n",
      "type", "cause", "benefit", "insurance_xj"
    )
  )

  v1 <- 1 / 1.04
  expected1 <- v1 * 0.10 + v1^2 * 0.70 * 0.20

  v2 <- 1 / 1.05
  expected2 <- 2 * (v2 * 0.20 + v2^2 * 0.70 * 0.30)

  expect_equal(
    out$insurance_xj,
    c(expected1, expected2),
    tolerance = 1e-12
  )
})

test_that("insurance_xj validates cause and term requirements", {
  md <- make_md_insurance_test()

  expect_error(
    insurance_xj(
      md = md,
      x = 30,
      i = 0.05,
      cause = "q_unknown",
      type = "term",
      n = 2
    ),
    "not found"
  )

  expect_error(
    insurance_xj(
      md = md,
      x = 30,
      i = 0.05,
      cause = "q_death",
      type = "term"
    ),
    "must be provided"
  )
})
