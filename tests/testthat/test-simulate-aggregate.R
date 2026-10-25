test_that("compound Poisson simulation reproduces aggregate mean", {
  model <- collective_risk_spec(
    frequency_spec("poisson", lambda = 2),
    severity_spec("degenerate", value = 100)
  )
  result <- simulate_aggregate(model, n_sim = 5000, seed = 100)
  expect_named(result, c("simulation", "claim_count", "aggregate_loss"))
  expect_true(all(result$aggregate_loss == 100 * result$claim_count))
  expect_true(all(result$aggregate_loss[result$claim_count == 0] == 0))
  expect_lt(abs(mean(result$aggregate_loss) - 200), 10)
})

test_that("compound binomial simulation aggregates independent severities", {
  model <- collective_risk_spec(
    frequency_spec("binomial", size = 2, prob = 1),
    severity_spec("discrete", values = c(10, 20), probs = c(0.5, 0.5))
  )
  out <- simulate_aggregate(model, 3000, seed = 14)
  expect_true(all(out$claim_count == 2))
  expect_true(all(out$aggregate_loss %in% c(20, 30, 40)))
  expect_lt(abs(mean(out$aggregate_loss) - 30), 0.7)
})

test_that("aggregate simulation supports exact zero-payment basis", {
  model <- collective_risk_spec(
    frequency_spec("poisson", lambda = 2),
    apply_coverage(severity_spec("degenerate", value = 100),
                   coverage_spec(coinsurance = 0))
  )
  zero <- to_payment_basis(model)
  expect_s3_class(zero, "tidyact_zero_payment_basis_spec")
  out <- simulate_aggregate(zero, 10, seed = 4)
  expect_true(all(out$claim_count == 0))
  expect_true(all(out$aggregate_loss == 0))
})

test_that("aggregate simulation restores seed and refuses silent truncation", {
  model <- collective_risk_spec(
    frequency_spec("binomial", size = 20, prob = 1),
    severity_spec("degenerate", value = 1)
  )
  set.seed(50)
  previous <- .Random.seed
  expect_identical(simulate_aggregate(model, 3, seed = 7),
                   simulate_aggregate(model, 3, seed = 7))
  expect_identical(.Random.seed, previous)
  expect_error(simulate_aggregate(model, n_sim = 3, max_claims = 50),
               "max_claims")
})

test_that("payment basis has equal analytic aggregate moments", {
  source <- severity_spec("discrete", values = c(0, 100), probs = c(0.5, 0.5))
  paid <- apply_coverage(source, coverage_spec(deductible = 0))
  per_loss <- collective_risk_spec(frequency_spec("poisson", lambda = 3), paid)
  per_payment <- to_payment_basis(per_loss)
  a <- aggregate_moments(per_loss, order = 1:2)
  b <- aggregate_moments(per_payment, order = 1:2)
  expect_equal(a$value, b$value)
  expect_equal(a$value[1], 150)
})

test_that("aggregate simulation can model both payment bases", {
  loss <- severity_spec("discrete", values = c(0, 100),
                        probs = c(0.25, 0.75))
  paid <- apply_coverage(loss, coverage_spec(deductible = 0))
  original <- collective_risk_spec(frequency_spec("poisson", lambda = 4), paid)
  payment_basis <- to_payment_basis(original)
  s1 <- simulate_aggregate(original, 4000, seed = 110)
  s2 <- simulate_aggregate(payment_basis, 4000, seed = 111)
  expect_lt(abs(mean(s1$aggregate_loss) - 300), 18)
  expect_lt(abs(mean(s2$aggregate_loss) - 300), 18)
})
