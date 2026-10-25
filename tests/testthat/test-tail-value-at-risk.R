test_that("TVaR allocates only the required probability at the VaR atom", {
  loss <- c(0, 100, 100, 300)
  result <- tail_value_at_risk(loss, level = c(0.5, 0.75, 0.9))
  expect_equal(result$tail_value_at_risk, c(200, 300, 300))
})

test_that("exponential TVaR equals VaR plus mean excess", {
  x <- severity_spec("exponential", rate = 0.01)
  p <- c(0.90, 0.95)
  expect_equal(tail_value_at_risk(x, p)$tail_value_at_risk,
               stats::qexp(p, rate = 0.01) + 100)
})

test_that("heavy-tailed severity with infinite mean has infinite TVaR", {
  x <- severity_spec("pareto", shape = 0.8, scale = 100)
  expect_equal(tail_value_at_risk(x, 0.95)$tail_value_at_risk, Inf)
})

test_that("TVaR refuses any positive unresolved upper-tail mass", {
  grid <- structure(list(grid = tibble::tibble(
    loss = c(0, 100), probability = c(0.5, 0.45)),
    tail_probability = 0.05), class = "tidyact_aggregate_distribution")
  expect_error(tail_value_at_risk(grid, 0.8), "unresolved upper tail")
})

test_that("TVaR handles capped continuous payments with a cap atom", {
  covered <- apply_coverage(
    severity_spec("exponential", rate = 0.01),
    coverage_spec(deductible = 50, limit = 200)
  )
  # At level .95, the upper cap has probability P(X >= 250) > .05.
  expect_equal(tail_value_at_risk(covered, level = 0.95)$tail_value_at_risk,
               200, tolerance = 1e-7)
})
