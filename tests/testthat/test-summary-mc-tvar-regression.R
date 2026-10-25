# Regression tests for coherent empirical VaR and TVaR.

test_that("TVaR allocates only the required fraction of an atom at VaR", {
  x <- c(0, 100, 100, 200)
  result <- summary_mc(tibble::tibble(loss = x), col_value = "loss",
                       probs = 0.5, var_probs = c(0.5, 0.75))
  expect_equal(result$VaR_500, 100)
  expect_equal(result$TVaR_500, 150)
  expect_equal(result$VaR_750, 100)
  expect_equal(result$TVaR_750, 200)
})

test_that("descriptive quantiles retain interpolation but empirical VaR does not", {
  x <- c(1, 2, 3, 4, 5)
  result <- summary_mc(tibble::tibble(loss = x), col_value = "loss",
                       probs = 0.8, var_probs = c(0.5, 0.8))
  expect_equal(result$q800, 4.2)
  expect_equal(result$VaR_800, 4)
  expect_equal(result$TVaR_500, 4.2)
  expect_equal(result$TVaR_800, 5)
})

test_that("empirical boundaries and ties have meaningful TVaR limits", {
  x <- c(0, 0, 100, 200)
  result <- summary_mc(tibble::tibble(loss = x), col_value = "loss",
                       var_probs = c(0, 0.5, 1))
  expect_equal(result$TVaR_000, mean(x))
  expect_equal(result$TVaR_500, 150)
  expect_equal(result$TVaR_1000, max(x))
  repeated <- summary_mc(tibble::tibble(loss = rep(40, 5)),
                         col_value = "loss", var_probs = c(0.1, 0.5, 0.9))
  expect_equal(unlist(repeated[1, c("TVaR_100", "TVaR_500", "TVaR_900")],
                      use.names = FALSE),
               c(40, 40, 40))
})

test_that("generic Monte Carlo summaries can retain negative values", {
  x <- c(-10, 0, 20, 50)
  result <- summary_mc(tibble::tibble(profit = x), col_value = "profit",
                       var_probs = 0.5)
  expect_equal(result$VaR_500, 0)
  expect_equal(result$TVaR_500, 35)
})

test_that("summary_mc agrees with standalone empirical actuarial measures", {
  x <- c(0, 100, 100, 200, 300)
  out <- summary_mc(tibble::tibble(loss = x), col_value = "loss",
                    var_probs = 0.8)
  expect_equal(out$VaR_800, value_at_risk(x, level = 0.8)$value_at_risk)
  expect_equal(out$TVaR_800, tail_value_at_risk(x, level = 0.8)$tail_value_at_risk,
               tolerance = 1e-12)
})
