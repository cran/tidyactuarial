.manual_bond_scenario_price <- function(face, c, n, k, j, R) {
  N <- as.integer(round(n * k))
  coupon <- face * c / k
  periods <- seq_len(N)
  sum(coupon / (1 + j)^periods) + R / (1 + j)^N
}

test_that("bond_book_value at time zero equals bond price", {
  j <- 0.03

  bv0 <- bond_book_value(
    face = 1000,
    c = 0.05,
    n = 10,
    t = 0,
    k = 2,
    y_effective_per_period = j
  )

  price <- .manual_bond_scenario_price(
    face = 1000,
    c = 0.05,
    n = 10,
    k = 2,
    j = j,
    R = 1000
  )

  expect_equal(bv0, price, tolerance = 1e-10)
})

test_that("par bond book value remains at face between coupon dates and is zero after maturity payment", {
  out <- bond_book_value(
    face = 1000,
    c = 0.06,
    n = 5,
    t = 0:5,
    k = 1,
    y_effective_per_period = 0.06
  )

  expect_equal(out[1:5], rep(1000, 5), tolerance = 1e-9)
  expect_equal(out[6], 0)
})

test_that("premium bond book value amortizes toward redemption", {
  out <- bond_book_value(
    face = 1000,
    c = 0.08,
    n = 5,
    t = 0:5,
    k = 1,
    y_effective_per_period = 0.05
  )

  expect_true(all(diff(out[1:5]) < 0))
  expect_true(all(out[1:5] > 1000))
  expect_equal(out[6], 0)
})

test_that("bond_book_value tidy output exposes valuation schedule", {
  out <- bond_book_value(
    face = 1000,
    c = 0.05,
    n = 4,
    t = c(0, 1, 2),
    k = 1,
    y_effective_per_period = 0.04,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$t, c(0, 1, 2))
  expect_equal(out$valuation_period, c(0L, 1L, 2L))
  expect_true(all(c("book_value", "yield_per_period", "yield_effective_annual") %in% names(out)))
})

test_that("bond_book_value requires valuation dates aligned with coupon dates", {
  expect_error(
    bond_book_value(
      face = 1000,
      c = 0.05,
      n = 5,
      t = 1.25,
      k = 2,
      y_effective_per_period = 0.03
    ),
    "valuation dates"
  )
})

test_that("bond_callable_price is the minimum price across call and maturity scenarios", {
  face <- 1000
  c <- 0.08
  n <- 10
  k <- 1
  j <- 0.05
  call_t <- c(5, 8)
  call_R <- c(1020, 1010)

  expected <- c(
    .manual_bond_scenario_price(face, c, call_t[1], k, j, call_R[1]),
    .manual_bond_scenario_price(face, c, call_t[2], k, j, call_R[2]),
    .manual_bond_scenario_price(face, c, n, k, j, face)
  )

  actual <- bond_callable_price(
    face = face,
    c = c,
    n = n,
    k = k,
    call_t = call_t,
    call_R = call_R,
    y_effective_per_period = j
  )

  expect_equal(actual, min(expected), tolerance = 1e-10)
})

test_that("bond_callable_price tidy output identifies the worst scenario", {
  out <- bond_callable_price(
    face = 1000,
    c = 0.08,
    n = 10,
    k = 1,
    call_t = c(5, 8),
    call_R = c(1020, 1010),
    y_effective_per_period = 0.05,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 3L)
  expect_equal(sum(out$is_worst_case), 1L)
  expect_equal(min(out$price_at_target_yield), bond_callable_price(
    face = 1000,
    c = 0.08,
    n = 10,
    k = 1,
    call_t = c(5, 8),
    call_R = c(1020, 1010),
    y_effective_per_period = 0.05
  ), tolerance = 1e-10)
  expect_equal(tail(out$scenario_type, 1), "maturity")
})

test_that("call dates must be strictly inside the original maturity", {
  expect_error(
    bond_callable_price(
      face = 1000,
      c = 0.05,
      n = 5,
      k = 1,
      call_t = 5,
      call_R = 1000,
      y_effective_per_period = 0.04
    ),
    "strictly between"
  )
})
