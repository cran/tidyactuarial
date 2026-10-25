test_that("annuity_arith increasing immediate matches direct discounted payments", {
  i <- 0.05
  n <- 6
  payments <- seq_len(n)
  times <- seq_len(n)
  expected <- sum(payments * (1 + i)^(-times))

  expect_equal(
    annuity_arith(n = n, i = i, pattern = "increasing", timing = "immediate"),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_arith decreasing due matches direct discounted payments", {
  i <- 0.04
  n <- 5
  payments <- n:1
  times <- 0:(n - 1)
  expected <- sum(payments * (1 + i)^(-times))

  expect_equal(
    annuity_arith(n = n, i = i, pattern = "decreasing", timing = "due"),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_arith custom pattern matches an explicit arithmetic stream", {
  i <- 0.06
  n <- 7
  P1 <- 100
  g <- 25
  payments <- P1 + (seq_len(n) - 1) * g
  expected <- sum(payments * (1 + i)^(-seq_len(n)))

  expect_equal(
    annuity_arith(
      n = n,
      i = i,
      pattern = "custom",
      P1 = P1,
      g = g
    ),
    expected,
    tolerance = 1e-10
  )
})

test_that("annuity_arith accumulated value matches direct accumulation", {
  i <- 0.05
  n <- 6
  payments <- seq_len(n)
  times <- seq_len(n)
  expected <- sum(payments * (1 + i)^(n - times))

  expect_equal(
    annuity_arith(
      n = n,
      i = i,
      pattern = "increasing",
      valuation = "accumulated"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_arith deferment affects present value but not terminal accumulated value", {
  i <- 0.05
  n <- 5
  h <- 2

  pv0 <- annuity_arith(n = n, i = i, h = 0, pattern = "increasing")
  pvh <- annuity_arith(n = n, i = i, h = h, pattern = "increasing")

  av0 <- annuity_arith(
    n = n,
    i = i,
    h = 0,
    pattern = "increasing",
    valuation = "accumulated"
  )
  avh <- annuity_arith(
    n = n,
    i = i,
    h = h,
    pattern = "increasing",
    valuation = "accumulated"
  )

  expect_equal(pvh, pv0 * (1 + i)^(-h), tolerance = 1e-12)
  expect_equal(avh, av0, tolerance = 1e-12)
})

test_that("annuity_arith handles fractional payment frequencies through the period rate", {
  i <- 0.1025
  k <- 2
  n <- 2
  i_period <- (1 + i)^(1 / k) - 1
  N <- n * k
  payments <- seq_len(N)
  expected <- sum(payments * (1 + i_period)^(-seq_len(N)))

  expect_equal(
    annuity_arith(n = n, k = k, i = i, pattern = "increasing"),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_arith tidy output exposes both valuation bases", {
  out <- annuity_arith(
    n = 5,
    i = 0.05,
    pattern = "custom",
    P1 = 50,
    g = 10,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_true(all(c(
    "present_value_factor",
    "accumulated_value_factor",
    "annuity_factor"
  ) %in% names(out)))
  expect_equal(out$annuity_factor, out$present_value_factor, tolerance = 1e-12)
})

test_that("annuity_geom finite immediate matches an explicit geometric stream", {
  i <- 0.06
  g <- 0.03
  n <- 8
  P1 <- 100
  payments <- P1 * (1 + g)^(0:(n - 1))
  expected <- sum(payments * (1 + i)^(-seq_len(n)))

  expect_equal(
    annuity_geom(n = n, i = i, g = g, P1 = P1),
    expected,
    tolerance = 1e-10
  )
})

test_that("annuity_geom due timing matches payments at times zero through n minus one", {
  i <- 0.05
  g <- 0.02
  n <- 7
  P1 <- 80
  payments <- P1 * (1 + g)^(0:(n - 1))
  expected <- sum(payments * (1 + i)^(-(0:(n - 1))))

  expect_equal(
    annuity_geom(n = n, i = i, g = g, P1 = P1, timing = "due"),
    expected,
    tolerance = 1e-10
  )
})

test_that("annuity_geom accumulated value matches direct accumulation", {
  i <- 0.07
  g <- 0.025
  n <- 6
  P1 <- 120
  payments <- P1 * (1 + g)^(0:(n - 1))
  times <- seq_len(n)
  expected <- sum(payments * (1 + i)^(n - times))

  expect_equal(
    annuity_geom(
      n = n,
      i = i,
      g = g,
      P1 = P1,
      valuation = "accumulated"
    ),
    expected,
    tolerance = 1e-10
  )
})

test_that("annuity_geom perpetuity matches the classical growing-perpetuity value", {
  i <- 0.08
  g <- 0.03
  P1 <- 250
  expected <- P1 / (i - g)

  expect_equal(
    annuity_geom(
      i = i,
      g = g,
      P1 = P1,
      perpetuity = TRUE,
      valuation = "present"
    ),
    expected,
    tolerance = 1e-10
  )
})

test_that("annuity_geom handles the equal interest-growth finite case", {
  i <- 0.05
  g <- 0.05
  n <- 9
  P1 <- 100
  expected <- n * P1 / (1 + i)

  expect_equal(
    annuity_geom(n = n, i = i, g = g, P1 = P1),
    expected,
    tolerance = 1e-10
  )
})

test_that("annuity_geom deferment discounts present value only", {
  i <- 0.06
  g <- 0.02
  n <- 6
  h <- 3

  pv0 <- annuity_geom(n = n, i = i, g = g, h = 0)
  pvh <- annuity_geom(n = n, i = i, g = g, h = h)

  av0 <- annuity_geom(n = n, i = i, g = g, h = 0, valuation = "accumulated")
  avh <- annuity_geom(n = n, i = i, g = g, h = h, valuation = "accumulated")

  expect_equal(pvh, pv0 * (1 + i)^(-h), tolerance = 1e-11)
  expect_equal(avh, av0, tolerance = 1e-11)
})

test_that("annuity_geom tidy output exposes present and accumulated values", {
  out <- annuity_geom(
    n = 5,
    i = 0.06,
    g = 0.02,
    P1 = 100,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_true(all(c(
    "present_value_factor",
    "accumulated_value_factor",
    "annuity_factor"
  ) %in% names(out)))
  expect_equal(out$annuity_factor, out$present_value_factor, tolerance = 1e-12)
})
