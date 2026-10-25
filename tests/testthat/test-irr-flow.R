test_that("irr_flow recovers a known one-period return", {
  result <- irr_flow(
    cf = c(-100, 110),
    t = c(0, 1),
    interval = c(0, 0.50)
  )

  expect_s3_class(result, "tbl_df")
  expect_true(result$converged)
  expect_equal(result$irr, 0.10, tolerance = 1e-9)
  expect_equal(result$npv, 0, tolerance = 1e-8)
  expect_equal(result$n_sign_changes, 1L)
})


test_that("irr_flow reports equivalent nominal interest and force", {
  result <- irr_flow(
    cf = c(-100, 121),
    t = c(0, 2),
    m = 12,
    interval = c(0, 0.50)
  )

  expected_irr <- 0.10
  expected_j12 <- 12 * ((1 + expected_irr)^(1 / 12) - 1)

  expect_equal(result$irr, expected_irr, tolerance = 1e-9)
  expect_equal(
    result$j_nominal_interest,
    expected_j12,
    tolerance = 1e-9
  )
  expect_equal(result$delta, log1p(expected_irr), tolerance = 1e-9)
})


test_that("irr_flow supports calendar dates", {
  result <- irr_flow(
    cf = c(-100, 110),
    date = as.Date(c("2026-01-01", "2027-01-01")),
    interval = c(0, 0.50)
  )

  expect_true(result$converged)
  expect_equal(result$irr, 0.10, tolerance = 1e-9)
})


test_that("irr_flow returns a diagnostic non-convergence result when no IRR can be bracketed", {
  same_sign <- irr_flow(
    cf = c(100, 50, 25),
    t = c(0, 1, 2)
  )

  expect_false(same_sign$converged)
  expect_true(is.na(same_sign$irr))
  expect_false(same_sign$has_both_signs)

  outside_interval <- irr_flow(
    cf = c(-100, 150),
    t = c(0, 1),
    interval = c(0, 0.20)
  )

  expect_false(outside_interval$converged)
  expect_true(is.na(outside_interval$irr))
  expect_true(outside_interval$has_both_signs)
})


test_that("irr_flow validates timing inputs", {
  expect_error(
    irr_flow(cf = c(-100, 110)),
    "either `t` or `date`"
  )

  expect_error(
    irr_flow(
      cf = c(-100, 110),
      t = c(0, 1),
      date = as.Date(c("2026-01-01", "2027-01-01"))
    ),
    "only one"
  )

  expect_error(
    irr_flow(
      cf = c(-100, 110),
      t = c(0, 1),
      interval = c(-1, 1)
    ),
    "greater than -1"
  )
})


test_that("irr_flow_multi detects two known actuarial roots", {
  # With x = 1 / (1 + r),
  # 50 - 115 x + 66 x^2 = 0
  # has roots x = 10/11 and x = 5/6,
  # corresponding to r = 0.10 and r = 0.20.
  result <- irr_flow_multi(
    cf = c(50, -115, 66),
    t = c(0, 1, 2),
    search_interval = c(0, 0.50),
    grid_points = 2001
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2L)
  expect_equal(result$irr, c(0.10, 0.20), tolerance = 1e-8)
  expect_equal(result$npv, c(0, 0), tolerance = 1e-7)
  expect_equal(result$n_sign_changes_cashflow, c(2L, 2L))
  expect_true(all(result$has_both_signs))
})


test_that("irr_flow_multi returns an empty tibble when a cash flow has one sign", {
  result <- irr_flow_multi(
    cf = c(100, 50, 25),
    t = c(0, 1, 2)
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0L)
  expect_named(
    result,
    c(
      "root_id",
      "irr",
      "i_effective_annual",
      "j_nominal_interest",
      "delta",
      "npv",
      "interval_left",
      "interval_right",
      "n_cashflows",
      "has_both_signs",
      "n_sign_changes_cashflow"
    )
  )
})


test_that("irr_flow_multi validates search controls", {
  expect_error(
    irr_flow_multi(
      cf = c(-100, 110),
      t = c(0, 1),
      search_interval = c(-1, 1)
    ),
    "greater than -1"
  )

  expect_error(
    irr_flow_multi(
      cf = c(-100, 110),
      t = c(0, 1),
      grid_points = 1
    ),
    "greater than or equal to 2"
  )
})
