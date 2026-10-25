test_that("Newton does not accept a remote point just because residual units are small", {
  for (scale in c(1, 1e-12, 1e-200, 1e200)) {
    out <- solve_value(fn = function(x) scale * (x - 2),
      solve_for = "x", target = 0, interval = c(0, 3), start = 0,
      derivative = function(x) scale, method = "newton", output = "audit")
    expect_equal(out$solution, 2, tolerance = 1e-10)
    expect_true(out$converged)
    estimated <- solve_value(fn = function(x) scale * (x - 2),
      solve_for = "x", target = 0, interval = c(0, 3), start = 0,
      method = "newton")
    expect_equal(estimated, 2, tolerance = 1e-10)
  }
})

test_that("a small residual at a constrained nonroot does not converge", {
  expect_error(solve_value(fn = function(x) 1e-12 * (x + 1),
    solve_for = "x", target = 0, interval = c(0, 2), start = 0,
    derivative = function(x) 1e-12, method = "newton"), "did not converge")
})
