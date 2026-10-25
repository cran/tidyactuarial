test_that("accrued_interest computes actual/actual coupon-period accrual", {
  x <- accrued_interest(
    face = 1000, c = 0.06, k = 2,
    previous_coupon = as.Date("2026-01-01"),
    settlement = as.Date("2026-04-01"),
    next_coupon = as.Date("2026-07-01"),
    clean_price = 1005
  )
  expected <- 30 * 90 / 181
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), 1L)
  expect_equal(x$coupon_per_period, 30)
  expect_equal(x$days_elapsed, 90)
  expect_equal(x$days_in_period, 181)
  expect_equal(x$year_fraction, 90 / (2 * 181), tolerance = 1e-12)
  expect_equal(x$coupon_fraction, 90 / 181, tolerance = 1e-12)
  expect_equal(x$accrued_interest, expected, tolerance = 1e-12)
  expect_equal(x$dirty_price, 1005 + expected, tolerance = 1e-12)
  expect_equal(x$clean_price, 1005)
})

test_that("accrued_interest accepts dirty price and computes clean price", {
  args <- list(face = 1000, c = 0.06, k = 2,
               previous_coupon = as.Date("2026-01-01"),
               settlement = as.Date("2026-04-01"),
               next_coupon = as.Date("2026-07-01"))
  x <- do.call(accrued_interest, c(args, list(dirty_price = 1020)))
  expect_equal(x$clean_price, 1020 - 30 * 90 / 181, tolerance = 1e-12)
  expect_equal(x$dirty_price, 1020)
})

test_that("accrued_interest applies act/360 and act/365 rather than prorating coupons", {
  args <- list(face = 1000, c = 0.06, k = 2,
               previous_coupon = as.Date("2026-01-01"),
               settlement = as.Date("2026-04-01"),
               next_coupon = as.Date("2026-07-01"))
  a360 <- do.call(accrued_interest, c(args, list(day_count = "act/360")))
  a365 <- do.call(accrued_interest, c(args, list(day_count = "act/365")))
  expect_equal(a360$accrued_interest, 1000 * 0.06 * 90 / 360)
  expect_equal(a365$accrued_interest, 1000 * 0.06 * 90 / 365)
  expect_gt(a360$accrued_interest, a365$accrued_interest)
  expect_equal(a360$coupon_fraction, 2 * 90 / 360)
})

test_that("accrued_interest implements European 30E/360 end-of-month cap", {
  x <- accrued_interest(
    face = 1000, c = 0.06, k = 2,
    previous_coupon = as.Date("2026-01-31"),
    settlement = as.Date("2026-04-30"),
    next_coupon = as.Date("2026-07-31"),
    day_count = "30e/360"
  )
  expect_equal(x$accrued_interest, 15)
  expect_equal(x$year_fraction, 0.25)
  expect_equal(x$coupon_fraction, 0.5)
})

test_that("accrued_interest correctly treats leap-year date arithmetic", {
  x <- accrued_interest(
    face = 100, c = 0.12, k = 1,
    previous_coupon = as.Date("2024-01-01"),
    settlement = as.Date("2024-03-01"),
    next_coupon = as.Date("2025-01-01")
  )
  expect_equal(x$days_in_period, 366)
  expect_equal(x$days_elapsed, 60)
  expect_equal(x$accrued_interest, 12 * 60 / 366)
})

test_that("accrued_interest resets to zero on previous coupon date", {
  dt <- as.Date("2026-01-01")
  x <- accrued_interest(face = 500, c = 0.04,
                        previous_coupon = dt, settlement = dt,
                        next_coupon = as.Date("2027-01-01"),
                        clean_price = 490)
  expect_equal(x$accrued_interest, 0)
  expect_equal(x$dirty_price, 490)
})

test_that("accrued_interest zero-coupon bond has zero accrued coupon interest", {
  x <- accrued_interest(
    face = 1000, c = 0, k = 2,
    previous_coupon = as.Date("2026-01-01"),
    settlement = as.Date("2026-05-01"),
    next_coupon = as.Date("2026-07-01")
  )
  expect_equal(x$accrued_interest, 0)
  expect_true(is.na(x$clean_price))
  expect_true(is.na(x$dirty_price))
})

test_that("accrued_interest validates both supplied prices consistently", {
  args <- list(face = 1000, c = 0.06, k = 2,
               previous_coupon = as.Date("2026-01-01"),
               settlement = as.Date("2026-04-01"),
               next_coupon = as.Date("2026-07-01"))
  amount <- 30 * 90 / 181
  x <- do.call(accrued_interest, c(args, list(
    clean_price = 1000, dirty_price = 1000 + amount)))
  expect_equal(x$dirty_price, 1000 + amount)
  expect_error(do.call(accrued_interest, c(args, list(
    clean_price = 1000, dirty_price = 1005))), "inconsistent")
})

test_that("accrued_interest works rowwise with a pipe-friendly bond tibble", {
  bonds <- tibble::tibble(
    bond_id = c("A", "B"), face = c(1000, 500), c = c(0.06, 0.04),
    k = c(2L, 2L),
    previous_coupon = as.Date(c("2026-01-01", "2026-01-01")),
    settlement = as.Date(c("2026-04-01", "2026-04-01")),
    next_coupon = as.Date(c("2026-07-01", "2026-07-01")),
    clean_price = c(1005, 500)
  )
  result <- bonds |> accrued_interest()
  expect_equal(result$bond_id, c("A", "B"))
  expect_equal(result$accrued_interest,
               c(30, 10) * 90 / 181, tolerance = 1e-12)
  expect_equal(result$dirty_price, result$clean_price + result$accrued_interest)
  expect_equal(nrow(result), 2L)
})

test_that("accrued_interest accepts per-bond day-count conventions", {
  data <- tibble::tibble(
    face = c(1000, 1000), c = c(0.06, 0.06),
    k = c(2L, 2L), day_count = c("act/act", "act/360"),
    previous_coupon = as.Date(c("2026-01-01", "2026-01-01")),
    settlement = as.Date(c("2026-04-01", "2026-04-01")),
    next_coupon = as.Date(c("2026-07-01", "2026-07-01"))
  )
  x <- accrued_interest(data)
  expect_equal(x$accrued_interest, c(30 * 90 / 181, 60 * 90 / 360))
})

test_that("accrued_interest supports vector inputs with scalar recycling", {
  x <- accrued_interest(face = c(1000, 500), c = 0.06, k = 2,
    previous_coupon = as.Date("2026-01-01"),
    settlement = as.Date("2026-04-01"),
    next_coupon = as.Date("2026-07-01"))
  expect_equal(x$accrued_interest, c(30, 15) * 90 / 181)
  expect_equal(nrow(x), 2L)
})

test_that("accrued_interest validates date order and coupon-date boundaries", {
  args <- list(face = 1000, c = 0.06, k = 2,
               previous_coupon = as.Date("2026-01-01"),
               settlement = as.Date("2026-04-01"),
               next_coupon = as.Date("2026-07-01"))
  expect_error(do.call(accrued_interest, modifyList(args, list(
    settlement = as.Date("2026-07-01")))), "Dates must satisfy")
  expect_error(do.call(accrued_interest, modifyList(args, list(
    settlement = as.Date("2025-12-31")))), "Dates must satisfy")
  expect_error(do.call(accrued_interest, modifyList(args, list(
    settlement = "2026-04-01"))), "valid Date")
})

test_that("accrued_interest rejects invalid financial or numerical inputs", {
  args <- list(face = 1000, c = 0.06, k = 2,
               previous_coupon = as.Date("2026-01-01"),
               settlement = as.Date("2026-04-01"),
               next_coupon = as.Date("2026-07-01"))
  expect_error(do.call(accrued_interest, modifyList(args, list(face = 0))),
               "strictly positive")
  expect_error(do.call(accrued_interest, modifyList(args, list(c = -0.01))),
               "nonnegative")
  expect_error(do.call(accrued_interest, modifyList(args, list(k = 2.5))),
               "positive integers")
  expect_error(do.call(accrued_interest, modifyList(args, list(day_count = "30/360"))),
               "day_count")
  expect_error(do.call(accrued_interest, modifyList(args, list(face = c(1000, 500, 200), c = c(0.06, 0.05)))),
               "length one or exactly")
  expect_error(do.call(accrued_interest, modifyList(args, list(tol = 0))),
               "positive numeric scalar")
})

test_that("accrued_interest rejects unexpected columns and missing input fields", {
  x <- tibble::tibble(face = 100, c = 0.05,
    settlement = as.Date("2026-02-01"),
    previous_coupon = as.Date("2026-01-01"),
    next_coupon = as.Date("2027-01-01"))
  expect_error(accrued_interest(x[, !names(x) %in% "c", drop = FALSE]),
               "Missing required column")
  expect_error(accrued_interest(dplyr::mutate(x, accrued_interest = 0)),
               "existing audit columns")
  expect_error(accrued_interest(x, face = 100), "read `face`")
  expect_error(accrued_interest(x, dirty_price = 100, clean_price = 95),
               "inconsistent")
})

test_that("accrued_interest coupon agrees with bond_cash_flows", {
  x <- accrued_interest(face = 1000, c = 0.06, k = 2,
    previous_coupon = as.Date("2026-01-01"),
    settlement = as.Date("2026-04-01"),
    next_coupon = as.Date("2026-07-01"))
  schedule <- bond_cash_flows(face = 1000, c = 0.06, n = 1, k = 2)
  expect_equal(x$coupon_per_period, schedule$cf[schedule$type == "coupon"][1])
})
