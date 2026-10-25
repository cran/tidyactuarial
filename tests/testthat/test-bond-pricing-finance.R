.manual_bond_price <- function(face, c, n, k = 1L, j, R = face) {
  N <- as.integer(round(n * k))
  if (N == 0L) return(R)
  coupon <- face * c / k
  periods <- seq_len(N)
  sum(coupon / (1 + j)^periods) + R / (1 + j)^N
}

test_that("bond_cash_flows builds the actuarial coupon and redemption schedule", {
  out <- bond_cash_flows(face = 1000, c = 0.06, n = 2, k = 2)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 5L)
  expect_equal(out$t, c(0.5, 1, 1.5, 2, 2))
  expect_equal(out$cf, c(30, 30, 30, 30, 1000))
  expect_equal(out$type, c("coupon", "coupon", "coupon", "coupon", "redemption"))
  expect_equal(sum(out$cf), 1120)
})

test_that("bond_cash_flows allows a redemption value different from face", {
  out <- bond_cash_flows(face = 1000, c = 0.05, n = 1, k = 1, R = 1020)
  expect_equal(out$cf, c(50, 1020))
  expect_equal(out$t, c(1, 1))
})

test_that("bond_cash_flows at maturity zero is immediate redemption only", {
  out <- bond_cash_flows(face = 1000, c = 0.05, n = 0, k = 2, R = 990)
  expect_equal(nrow(out), 1L)
  expect_equal(out$t, 0)
  expect_equal(out$cf, 990)
  expect_equal(out$type, "redemption")
})

test_that("bond_cash_flows rejects stub periods", {
  expect_error(
    bond_cash_flows(face = 1000, c = 0.05, n = 1.25, k = 2),
    "n \\* k"
  )
})

test_that("bond_price agrees with direct present-value calculation", {
  face <- 1000
  c <- 0.07
  n <- 6
  k <- 2
  j <- 0.025
  R <- 1010

  expected <- .manual_bond_price(face, c, n, k, j, R)
  actual <- bond_price(
    face = face,
    c = c,
    n = n,
    k = k,
    y_effective_per_period = j,
    R = R
  )

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("bond_price equals par when coupon rate matches nominal yield frequency", {
  expect_equal(
    bond_price(
      face = 1000,
      c = 0.06,
      n = 10,
      k = 2,
      y_effective_per_period = 0.03
    ),
    1000,
    tolerance = 1e-9
  )
})

test_that("bond_price annual-effective yield agrees with equivalent per-period yield", {
  j <- 0.03
  i_annual <- (1 + j)^2 - 1

  direct <- bond_price(
    face = 1000,
    c = 0.05,
    n = 8,
    k = 2,
    y_effective_per_period = j
  )

  annual <- bond_price(
    face = 1000,
    c = 0.05,
    n = 8,
    k = 2,
    y = i_annual,
    y_type = "effective"
  )

  expect_equal(annual, direct, tolerance = 1e-10)
})

test_that("bond_price responds correctly to yield relative to coupon rate", {
  premium <- bond_price(
    face = 1000, c = 0.08, n = 5, k = 1,
    y_effective_per_period = 0.05
  )
  discount <- bond_price(
    face = 1000, c = 0.04, n = 5, k = 1,
    y_effective_per_period = 0.06
  )

  expect_gt(premium, 1000)
  expect_lt(discount, 1000)
})

test_that("bond_price at zero remaining term equals redemption value", {
  expect_equal(
    bond_price(
      face = 1000,
      c = 0.06,
      n = 0,
      k = 2,
      y_effective_per_period = 0.03,
      R = 1025
    ),
    1025
  )
})

test_that("bond_ytm inverts bond_price", {
  j <- 0.0275
  P <- .manual_bond_price(
    face = 1000,
    c = 0.055,
    n = 9,
    k = 2,
    j = j,
    R = 1000
  )

  out <- bond_ytm(
    P = P,
    face = 1000,
    c = 0.055,
    n = 9,
    k = 2
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$i_period, j, tolerance = 1e-9)
  expect_equal(out$j_nominal, 2 * j, tolerance = 1e-9)
  expect_equal(out$i_effective_annual, (1 + j)^2 - 1, tolerance = 1e-9)
})

test_that("bond_ytm recovers coupon yield for a par annual bond", {
  out <- bond_ytm(
    P = 1000,
    face = 1000,
    c = 0.07,
    n = 8,
    k = 1
  )

  expect_equal(out$i_period, 0.07, tolerance = 1e-9)
  expect_equal(out$i_effective_annual, 0.07, tolerance = 1e-9)
})

test_that("bond_ytm validates supplied root interval", {
  expect_error(
    bond_ytm(
      P = 1000,
      face = 1000,
      c = 0.05,
      n = 5,
      interval = c(0.10, 0.05)
    ),
    "interval"
  )
})
