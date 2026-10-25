test_that("frequency simulation returns nonnegative counts with correct moments", {
  x <- frequency_spec("poisson", lambda = 3)
  result <- simulate_frequency(x, n_sim = 6000, seed = 108)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("simulation", "claim_count"))
  expect_equal(result$simulation, seq_len(6000))
  expect_true(all(result$claim_count >= 0))
  expect_true(all(result$claim_count == floor(result$claim_count)))
  expect_lt(abs(mean(result$claim_count) - 3), 0.2)
})

test_that("frequency simulation supports degenerate binomial and empirical frequency", {
  expect_true(all(simulate_frequency(
    frequency_spec("binomial", size = 30, prob = 0), 50, seed = 2
  )$claim_count == 0))
  empirical <- empirical_frequency(c(0, 2, 2, 4))
  result <- simulate_frequency(empirical, 100, seed = 21)
  expect_true(all(result$claim_count %in% c(0, 2, 4)))
})

test_that("frequency simulation preserves caller random state when seeded", {
  set.seed(17)
  previous <- .Random.seed
  x <- frequency_spec("negative_binomial", size = 2, prob = 0.3)
  first <- simulate_frequency(x, 20, seed = 101)
  expect_identical(.Random.seed, previous)
  expect_identical(first, simulate_frequency(x, 20, seed = 101))
})

test_that("frequency simulation rejects unavailable RNGs and invalid arguments", {
  custom <- frequency_spec("custom", custom = list(
    pmf = function(x, params, log = FALSE) stats::dpois(x, lambda = 2, log = log)
  ))
  expect_error(simulate_frequency(custom, 2), "no RNG provider")
  expect_error(simulate_frequency(frequency_spec("poisson", lambda = 1), 0), "n_sim")
  expect_error(simulate_frequency(frequency_spec("poisson", lambda = 1), 2,
                                  seed = -1), "seed")
})
