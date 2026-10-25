test_that("standardize_interest implements the four standard conversions", {
  j12 <- 0.12
  d4 <- 0.08
  delta <- 0.05

  expect_equal(
    standardize_interest(i_type = "effective", i = 0.07),
    0.07
  )

  expect_equal(
    standardize_interest(
      i_type = "nominal_interest",
      i = j12,
      m = 12
    ),
    (1 + j12 / 12)^12 - 1
  )

  expect_equal(
    standardize_interest(
      i_type = "nominal_discount",
      i = d4,
      m = 4
    ),
    (1 - d4 / 4)^(-4) - 1
  )

  expect_equal(
    standardize_interest(i_type = "force", i = delta),
    exp(delta) - 1
  )
})


test_that("standardize_interest vectorizes heterogeneous rate conventions", {
  i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )
  i <- c(0.06, 0.12, 0.08, 0.05)
  m <- c(1, 12, 4, 1)

  expected <- c(
    0.06,
    (1 + 0.12 / 12)^12 - 1,
    (1 - 0.08 / 4)^(-4) - 1,
    exp(0.05) - 1
  )

  expect_equal(
    standardize_interest(i_type = i_type, i = i, m = m),
    expected
  )
})


test_that("standardize_interest preserves transitional type/rate compatibility", {
  expect_equal(
    standardize_interest(
      type = "nominal_interest",
      rate = 0.12,
      m = 12
    ),
    (1 + 0.12 / 12)^12 - 1
  )
})


test_that("standardize_interest propagates missing rate values", {
  expect_equal(
    standardize_interest(
      i_type = c("effective", "force", "effective"),
      i = c(0.05, NA, 0.07),
      m = 1
    ),
    c(0.05, NA_real_, 0.07)
  )
})


test_that("standardize_interest rejects invalid rate specifications", {
  expect_error(
    standardize_interest(i_type = "unknown", i = 0.05),
    "must contain only"
  )

  expect_error(
    standardize_interest(
      i_type = "nominal_interest",
      i = 0.12,
      m = 0
    ),
    "positive integer"
  )

  expect_error(
    standardize_interest(
      i_type = "effective",
      i = -1
    ),
    "greater than -1"
  )

  expect_error(
    standardize_interest(
      i_type = "nominal_discount",
      i = 4,
      m = 4
    ),
    "must be positive"
  )

  expect_error(
    standardize_interest(
      i = 0.05,
      rate = 0.06
    ),
    "only one"
  )
})


test_that("interest_equivalents reproduces the classical equivalence identities", {
  i_eff <- 0.08
  m <- 12

  result <- interest_equivalents(
    i_type = "effective",
    i = i_eff,
    m = m
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 6L)
  expect_equal(
    result$family,
    c(
      "effective",
      "discount",
      "discount_factor",
      "force",
      "nominal_interest",
      "nominal_discount"
    )
  )

  expected <- c(
    i_eff,
    i_eff / (1 + i_eff),
    1 / (1 + i_eff),
    log1p(i_eff),
    m * ((1 + i_eff)^(1 / m) - 1),
    m * (1 - (1 + i_eff)^(-1 / m))
  )

  expect_equal(result$value, expected, tolerance = 1e-12)
})


test_that("interest_equivalents is invariant to the equivalent input convention", {
  i_eff <- 0.08
  delta <- log1p(i_eff)
  j12 <- 12 * ((1 + i_eff)^(1 / 12) - 1)

  from_effective <- interest_equivalents(
    i_type = "effective",
    i = i_eff,
    m = 12
  )

  from_force <- interest_equivalents(
    i_type = "force",
    i = delta,
    m = 12
  )

  from_nominal <- interest_equivalents(
    i_type = "nominal_interest",
    i = j12,
    m = 12
  )

  expect_equal(from_force$value, from_effective$value, tolerance = 1e-12)
  expect_equal(from_nominal$value, from_effective$value, tolerance = 1e-12)
})


test_that("interest_equivalents validates scalar inputs", {
  expect_error(
    interest_equivalents(i = c(0.05, 0.06)),
    "single finite numeric"
  )

  expect_error(
    interest_equivalents(i = 0.05, m = 2.5),
    "positive integer"
  )
})
