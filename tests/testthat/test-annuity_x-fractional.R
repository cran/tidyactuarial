test_that("annuity_x exact k-thly valuation uses payments of 1 / k", {
  lt <- data.frame(
    x = 60:66,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
  )

  i <- 0.05
  k <- 12L
  n <- 2

  lx_at <- function(age) {
    if (age == 67) {
      return(0)
    }

    lt$lx[match(age, lt$x)]
  }

  survival_udd <- function(u) {
    completed_years <- floor(u)
    fractional_year <- u - completed_years

    survival_to_year <- lx_at(60 + completed_years) / lx_at(60)

    if (abs(fractional_year) < 1e-12) {
      return(survival_to_year)
    }

    p_year <- lx_at(61 + completed_years) / lx_at(60 + completed_years)
    q_year <- 1 - p_year

    survival_to_year * (1 - fractional_year * q_year)
  }

  times <- (0:(k * n - 1L)) / k

  expected <- sum(
    (1 / k) *
      (1 + i)^(-times) *
      vapply(times, survival_udd, numeric(1L))
  )

  observed <- annuity_x(
    lt = lt,
    x = 60,
    i = i,
    n = n,
    k = k,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("annuity_x accepts fractional exact terms aligned with k", {
  lt <- data.frame(
    x = 60:66,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
  )

  result <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 2.5,
    k = 12,
    timing = "due",
    woolhouse = "none",
    frac = "UDD",
    tidy = TRUE
  )

  expect_equal(result$n_used, 2.5)
  expect_true(is.finite(result$apv))
  expect_gt(result$apv, 0)
})


test_that("annuity_x rejects exact terms not aligned with k", {
  lt <- data.frame(
    x = 60:66,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
  )

  expect_error(
    annuity_x(
      lt = lt,
      x = 60,
      i = 0.05,
      n = 2.3,
      k = 4,
      timing = "due",
      woolhouse = "none"
    ),
    "`n \\* k` must be an integer"
  )
})


test_that("annuity_x keeps Woolhouse terms integer-valued", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  expect_error(
    annuity_x(
      lt = lt,
      x = 60,
      i = 0.05,
      n = 2.5,
      k = 12,
      timing = "due",
      woolhouse = "first"
    ),
    "Woolhouse approximations require `n` to be an integer"
  )
})


test_that("annuity_x preserves the k-thly due-immediate identity", {
  lt <- data.frame(
    x = 60:66,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
  )

  i <- 0.05
  n <- 2.5
  k <- 12

  due <- annuity_x(
    lt = lt,
    x = 60,
    i = i,
    n = n,
    k = k,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  immediate <- annuity_x(
    lt = lt,
    x = 60,
    i = i,
    n = n,
    k = k,
    timing = "immediate",
    woolhouse = "none",
    frac = "UDD"
  )

  lx_60 <- lt$lx[lt$x == 60]
  lx_62 <- lt$lx[lt$x == 62]
  lx_63 <- lt$lx[lt$x == 63]

  p_2 <- lx_62 / lx_60
  p_62 <- lx_63 / lx_62
  q_62 <- 1 - p_62
  survival_2_5 <- p_2 * (1 - 0.5 * q_62)

  pure_endowment <- (1 + i)^(-n) * survival_2_5

  expect_equal(
    immediate,
    due - (1 / k) * (1 - pure_endowment),
    tolerance = 1e-12
  )
})
