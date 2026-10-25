test_that("coordinate accuracy is preserved when residual units change", {
  for (scale in c(1, 1e-12, 1e12)) {
    for (method in c("uniroot", "newton")) {
      for (explicit in c(FALSE, TRUE)) {
        out <- solve_value(fn = function(x) scale * x^2,
          solve_for = "x", target = 2 * scale, interval = c(1, 2),
          start = 1.2, derivative = if (explicit) function(x) 2 * scale * x else NULL,
          method = method, scan_points = 3L, tol = 1e-11,
          residual_tol = scale * 1e-9, output = "audit")
        expect_equal(out$solution, sqrt(2), tolerance = 1e-10)
        expect_lte(abs(out$achieved - 2 * scale), scale * 1e-9)
        expect_true(out$converged)
      }
    }
  }
})

test_that("audit uses the explicit residual threshold without a hidden floor", {
  args <- list(fn = function(x) 1e5 * (x^2 - 2), solve_for = "x",
    target = 0, interval = c(1, 2), scan_points = 3L, method = "uniroot",
    tol = 1e-14, output = "audit")
  implicit <- do.call(solve_value, args)
  strict <- do.call(solve_value, c(args, list(residual_tol = 1e-12)))
  loose <- do.call(solve_value, c(args, list(residual_tol = 1e-8)))
  expect_gt(abs(strict$residual), 1e-12)
  expect_lt(abs(strict$residual), 1e-8)
  expect_false(implicit$converged)
  expect_false(strict$converged)
  expect_true(loose$converged)
  expect_identical(strict$solution, loose$solution)
  # Retain the documented candidate return for value and summary outputs.
  args$output <- "value"
  expect_identical(do.call(solve_value, args), implicit$solution)
  args$output <- "summary"
  expect_identical(do.call(solve_value, args)$solution, implicit$solution)
})

test_that("existing positional arguments and default residual tolerance agree", {
  fn <- function(x) x^2 - 2
  positional <- solve_value(fn, "x", 0, list(), NULL, c(1, 2), NULL, NULL,
    "uniroot", "error", 1e-10, 1000L, 3L, 0L, "audit")
  named <- solve_value(fn, solve_for = "x", target = 0, interval = c(1, 2),
    method = "uniroot", tol = 1e-10, maxiter = 1000L, scan_points = 3L,
    max_expand = 0L, output = "audit", residual_tol = 1e-10)
  expect_identical(positional, named)
  zero <- solve_value(function(x) x, "x", 0, interval = c(-1, 1),
    residual_tol = 1e-20, output = "audit")
  expect_identical(zero$solution, 0)
  expect_true(zero$converged)
})

test_that("residual tolerance is validated and iteration exhaustion is an error", {
  for (bad in list(NULL, NA_real_, Inf, 0, -1, c(1, 2), "small", TRUE)) {
    expect_error(solve_value(function(x) x, "x", 0, residual_tol = bad),
      "`residual_tol` must be a single finite positive number", fixed = TRUE)
  }
  expect_error(solve_value(function(x) x^3 - 2, "x", 0,
    interval = c(0, 2), method = "uniroot", maxiter = 1L, scan_points = 3L,
    tol = 1e-14, output = "audit"), "convergence|iteration")
})
