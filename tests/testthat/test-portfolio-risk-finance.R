test_that("portfolio_duration is market-value weighted duration", {
  P <- c(100, 200, 300)
  D <- c(2, 5, 8)
  expected <- sum(P * D) / sum(P)

  out <- portfolio_duration(P = P, D = D)

  expect_s3_class(out, "tbl_df")
  expect_equal(out$D_P, expected, tolerance = 1e-12)
  expect_equal(out$P_total, 600)
  expect_equal(out$n_positions, 3L)
})

test_that("portfolio_convexity is market-value weighted convexity", {
  P <- c(100, 200, 300)
  C <- c(6, 30, 72)
  expected <- sum(P * C) / sum(P)

  out <- portfolio_convexity(P = P, C = C)

  expect_s3_class(out, "tbl_df")
  expect_equal(out$C_P, expected, tolerance = 1e-12)
  expect_equal(out$P_total, 600)
  expect_equal(out$n_positions, 3L)
})

test_that("portfolio duration and convexity support grouped tidy data", {
  dat <- tibble::tibble(
    book = c("A", "A", "B", "B"),
    value = c(100, 300, 200, 200),
    dur = c(2, 6, 4, 8),
    conv = c(5, 40, 20, 70)
  )

  d_out <- portfolio_duration(
    .data = dat,
    col_portfolio = "book",
    col_P = "value",
    col_D = "dur"
  )
  c_out <- portfolio_convexity(
    .data = dat,
    col_portfolio = "book",
    col_P = "value",
    col_C = "conv"
  )

  expect_equal(d_out$book, c("A", "B"))
  expect_equal(d_out$D_P, c(5, 6), tolerance = 1e-12)
  expect_equal(c_out$book, c("A", "B"))
  expect_equal(c_out$C_P, c(31.25, 45), tolerance = 1e-12)
})

test_that("portfolio risk functions allow custom output names", {
  d_out <- portfolio_duration(
    P = c(100, 200),
    D = c(3, 7),
    .out = "duration",
    .out_value = "market_value",
    .out_n = "positions"
  )

  c_out <- portfolio_convexity(
    P = c(100, 200),
    C = c(10, 40),
    .out = "convexity",
    .out_value = "market_value",
    .out_n = "positions"
  )

  expect_true(all(c("duration", "market_value", "positions") %in% names(d_out)))
  expect_true(all(c("convexity", "market_value", "positions") %in% names(c_out)))
})

test_that("portfolio risk functions propagate missing values by default", {
  d_out <- portfolio_duration(P = c(100, NA), D = c(3, 7))
  c_out <- portfolio_convexity(P = c(100, NA), C = c(10, 40))

  expect_true(is.na(d_out$D_P))
  expect_true(is.na(c_out$C_P))
})

test_that("portfolio risk functions reject negative market values", {
  expect_error(
    portfolio_duration(P = c(100, -20), D = c(3, 7)),
    "must be finite and >= 0"
  )
  expect_error(
    portfolio_convexity(P = c(100, -20), C = c(10, 40)),
    "must be finite and >= 0"
  )
})
