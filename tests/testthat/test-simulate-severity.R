test_that("severity simulation returns iid nonnegative claim amounts", {
  x <- severity_spec("exponential", rate = 0.01)
  result <- simulate_severity(x, n_sim = 5000, seed = 4)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("simulation", "claim_amount"))
  expect_true(all(result$claim_amount >= 0))
  expect_lt(abs(mean(result$claim_amount) - 100), 5)
})

test_that("severity simulation retains zeros from discrete and covered losses", {
  x <- severity_spec("discrete", values = c(0, 100), probs = c(0.8, 0.2))
  expect_true(any(simulate_severity(x, 100, seed = 100)$claim_amount == 0))
  covered <- apply_coverage(x, coverage_spec(deductible = 100))
  expect_true(all(simulate_severity(covered, 20, seed = 2)$claim_amount == 0))
  empirical <- empirical_severity(c(0, 50, 50))
  expect_true(all(simulate_severity(empirical, 30, seed = 8)$claim_amount %in%
                  c(0, 50)))
})

test_that("severity simulation restores RNG and validates providers", {
  set.seed(10)
  previous <- .Random.seed
  x <- severity_spec("gamma", shape = 2, rate = 0.02)
  expect_identical(simulate_severity(x, 10, seed = 9),
                   simulate_severity(x, 10, seed = 9))
  expect_identical(.Random.seed, previous)
  custom <- severity_spec("custom", custom = list(
    pdf = function(x, params, log = FALSE) stats::dexp(x, rate = 1, log = log),
    support = list(type = "continuous", lower = 0, upper = Inf)
  ))
  expect_error(simulate_severity(custom, 3), "no RNG provider")
  expect_error(simulate_severity(x, n_sim = -2), "n_sim")
})
