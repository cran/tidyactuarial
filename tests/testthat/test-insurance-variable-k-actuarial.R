make_lt_variable_insurance <- function() {
  lifetable(
    x = 60:64,
    lx = c(1000, 800, 600, 300, 0),
    close = TRUE,
    frac = "UDD"
  )
}

test_that("insurance_variable_k with k=1 matches an explicit annual death-benefit sum", {
  lt <- make_lt_variable_insurance()
  i <- 0.05
  v <- 1 / (1 + i)

  # Survival probabilities: 1, .8, .6, .3.
  # Death probabilities over the first 3 years: .2, .2, .3.
  expected <- v * 0.20 + v^2 * 0.20 + v^3 * 0.30

  expect_equal(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = i,
      benefit = 1,
      n = 3,
      k = 1
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("insurance_variable_k matches a manual UDD half-year calculation", {
  lt <- make_lt_variable_insurance()
  i <- 0.05
  v <- 1 / (1 + i)

  # Under UDD:
  # S(0)=1, S(.5)=.9, S(1)=.8, S(1.5)=.7, S(2)=.6.
  # Hence each half-year death probability here equals .1.
  expected <- 0.10 * (
    v^0.5 + v^1 + v^1.5 + v^2
  )

  expect_equal(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = i,
      benefit = 1,
      n = 2,
      k = 2,
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("insurance_variable_k handles genuinely varying subperiod benefits", {
  lt <- make_lt_variable_insurance()
  i <- 0.04
  v <- 1 / (1 + i)
  b <- c(100, 200, 300, 400)

  expected <- 0.10 * (
    100 * v^0.5 +
    200 * v^1 +
    300 * v^1.5 +
    400 * v^2
  )

  expect_equal(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = i,
      benefit = b,
      n = 2,
      k = 2,
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-10
  )
})

test_that("insurance_variable_k accepts a benefit function evaluated at payment times", {
  lt <- make_lt_variable_insurance()
  i <- 0.04
  v <- 1 / (1 + i)

  benefit_fun <- function(t) 1000 * t

  # At t=.5, 1, 1.5, 2; each death probability is .1 under UDD.
  expected <- 0.10 * (
    500 * v^0.5 +
    1000 * v^1 +
    1500 * v^1.5 +
    2000 * v^2
  )

  expect_equal(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = i,
      benefit = benefit_fun,
      n = 2,
      k = 2,
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-10
  )
})

test_that("insurance_variable_k applies deferment before the insured term", {
  lt <- make_lt_variable_insurance()
  i <- 0.05
  v <- 1 / (1 + i)

  # h=1, n=2, k=1 -> death during years 2 and 3.
  expected <- v^2 * (0.80 - 0.60) +
    v^3 * (0.60 - 0.30)

  expect_equal(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = i,
      benefit = 1,
      n = 2,
      h = 1,
      k = 1
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("insurance_variable_k infers n from a numeric benefit vector", {
  lt <- make_lt_variable_insurance()
  b <- c(100, 200, 300, 400)

  inferred <- insurance_variable_k(
    lt = lt,
    x = 60,
    i = 0.04,
    benefit = b,
    k = 2
  )

  explicit <- insurance_variable_k(
    lt = lt,
    x = 60,
    i = 0.04,
    benefit = b,
    n = 2,
    k = 2
  )

  expect_equal(inferred, explicit, tolerance = 1e-12)
})

test_that("insurance_variable_k respects equivalent interest conventions", {
  lt <- make_lt_variable_insurance()

  i_nominal <- 0.06
  m <- 12
  i_eff <- (1 + i_nominal / m)^m - 1

  nominal <- insurance_variable_k(
    lt = lt,
    x = 60,
    i = i_nominal,
    i_type = "nominal_interest",
    m = m,
    benefit = 1,
    n = 2,
    k = 2
  )

  effective <- insurance_variable_k(
    lt = lt,
    x = 60,
    i = i_eff,
    i_type = "effective",
    benefit = 1,
    n = 2,
    k = 2
  )

  expect_equal(nominal, effective, tolerance = 1e-12)
})

test_that("insurance_variable_k tidy output records the actuarial assumptions", {
  lt <- make_lt_variable_insurance()

  out <- insurance_variable_k(
    lt = lt,
    x = 60,
    i = 0.05,
    benefit = rep(1000, 4),
    n = 2,
    k = 2,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
  expect_named(
    out,
    c(
      "x", "h", "n", "k", "i", "i_type",
      "m", "i_effective", "frac", "apv"
    )
  )
  expect_equal(out$frac, "UDD")
  expect_true(is.finite(out$apv))
  expect_gt(out$apv, 0)
})

test_that("insurance_variable_k validates subperiod grids and benefit length", {
  lt <- make_lt_variable_insurance()

  expect_error(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = 0.05,
      benefit = 1,
      n = 1.1,
      k = 2
    ),
    "n \\* k"
  )

  expect_error(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = 0.05,
      benefit = c(100, 200, 300),
      n = 2,
      k = 2
    ),
    "length"
  )

  expect_error(
    insurance_variable_k(
      lt = lt,
      x = 60,
      i = 0.05,
      benefit = function(t) t,
      k = 2
    ),
    "Provide `n`"
  )
})
