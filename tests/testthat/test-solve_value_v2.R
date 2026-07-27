test_that("solve_value solves an effective interest rate", {
  expected <- (10000 / 8000)^(1 / 5) - 1

  result <- solve_value(
    fn = present_value,
    solve_for = "i",
    target = 8000,
    args = list(
      C = 10000,
      t = 5,
      i_type = "effective"
    ),
    interval = c(0, 0.20)
  )

  expect_equal(result, expected, tolerance = 1e-8)
})


test_that("solve_value extracts a result from a tidy output", {
  factor <- a_angle(
    n = 10,
    i = 0.05,
    timing = "immediate"
  )
  expected_payment <- 100000 / factor

  result <- solve_value(
    fn = a_angle,
    solve_for = "payment",
    target = 100000,
    args = list(
      n = 10,
      i = 0.05,
      timing = "immediate",
      tidy = TRUE
    ),
    result = "present_value",
    interval = c(0, 20000)
  )

  expect_equal(result, expected_payment, tolerance = 1e-8)
})


test_that("solve_value solves for a principal amount", {
  result <- solve_value(
    fn = present_value,
    solve_for = "C",
    target = 8000,
    args = list(
      i = 0.05,
      t = 5,
      i_type = "effective"
    ),
    interval = c(0, 20000)
  )

  expect_equal(result, 8000 * 1.05^5, tolerance = 1e-7)
})


test_that("solve_value can return all bracketed roots", {
  result <- solve_value(
    fn = function(x) (x - 1) * (x - 3),
    solve_for = "x",
    target = 0,
    interval = c(0, 4),
    multiple = "all",
    scan_points = 501
  )

  expect_equal(result, c(1, 3), tolerance = 1e-7)
})


test_that("solve_value detects multiple roots by default", {
  expect_error(
    solve_value(
      fn = function(x) (x - 1) * (x - 3),
      solve_for = "x",
      target = 0,
      interval = c(0, 4),
      scan_points = 501
    ),
    "Multiple roots"
  )
})


test_that("solve_value supports Newton-Raphson", {
  result <- solve_value(
    fn = function(x) x^2,
    solve_for = "x",
    target = 2,
    interval = c(0, 2),
    start = 1,
    method = "newton"
  )

  expect_equal(result, sqrt(2), tolerance = 1e-8)
})


test_that("summary output is compact", {
  result <- solve_value(
    fn = present_value,
    solve_for = "i",
    target = 8000,
    args = list(C = 10000, t = 5),
    interval = c(0, 0.20),
    output = "summary"
  )

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c("solve_for", "solution", "target", "achieved", "residual", "method")
  )
  expect_equal(ncol(result), 6L)
})


test_that("audit output reports diagnostics", {
  result <- solve_value(
    fn = present_value,
    solve_for = "i",
    target = 8000,
    args = list(C = 10000, t = 5),
    interval = c(0, 0.20),
    output = "audit"
  )

  expect_s3_class(result, "tbl_df")
  expect_true(all(c(
    "solution",
    "residual",
    "converged",
    "interval_left",
    "interval_right",
    "n_roots"
  ) %in% names(result)))
  expect_true(result$converged)
})


test_that("solve_value requires an extractor for non-scalar outputs", {
  expect_error(
    solve_value(
      fn = a_angle,
      solve_for = "payment",
      target = 100000,
      args = list(
        n = 10,
        i = 0.05,
        tidy = TRUE
      ),
      interval = c(0, 20000)
    ),
    "Use `result`"
  )
})


test_that("solve_value fails clearly when the target is unattainable", {
  expect_error(
    solve_value(
      fn = function(x) x^2 + 1,
      solve_for = "x",
      target = 0,
      interval = c(-2, 2),
      max_expand = 0
    ),
    "No root was found"
  )
})


test_that("an explicit interval is not expanded", {
  expect_error(
    solve_value(
      fn = function(x) x,
      solve_for = "x",
      target = 10,
      interval = c(0, 1)
    ),
    "No root was found"
  )
})


test_that("audit metadata stays aligned with multiple roots", {
  result <- solve_value(
    fn = function(x) (x - 1) * (x - 3),
    solve_for = "x",
    target = 0,
    interval = c(0, 4),
    multiple = "all",
    scan_points = 501,
    output = "audit"
  )

  expect_equal(result$solution, c(1, 3), tolerance = 1e-7)
  expect_equal(length(result$method), 2L)
  expect_equal(length(result$n_iter), 2L)
})
