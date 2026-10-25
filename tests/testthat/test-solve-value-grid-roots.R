test_that("small residual units do not turn the solve grid into many roots", {
  for (scale in c(1, 1e-12, 1e-200)) {
    out <- solve_value(fn = function(x) scale * (x - 1.137),
      solve_for = "x", target = 0, interval = c(0, 2),
      multiple = "all", method = "uniroot")
    expect_length(out, 1L)
    expect_equal(out, 1.137, tolerance = 1e-9)
  }
})

test_that("an unattainable small target residual is not a grid root", {
  expect_error(solve_value(fn = function(x) 1e-12 * (x + 1),
    solve_for = "x", target = 0, interval = c(0, 2),
    method = "uniroot"), "No root was found")
})
