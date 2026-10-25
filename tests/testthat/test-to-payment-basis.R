test_that("Poisson thinning and aggregate moments are invariant", {
  x <- severity_spec("discrete", values = c(50, 150, 250),
                     probs = c(.2, .5, .3))
  payment <- apply_coverage(x, coverage_spec(deductible = 100))
  original <- collective_risk_spec(frequency_spec("poisson", lambda = 4),
                                   payment)
  paid <- to_payment_basis(original)
  expect_s3_class(paid, "tidyact_collective_risk_spec")
  expect_identical(paid$count_basis, "payment")
  expect_equal(paid$thinning$probability, .8)
  expect_equal(paid$frequency$params$lambda, 3.2)
  expect_equal(paid$severity$params$values, c(50, 150))
  expect_equal(paid$severity$params$probs, c(.625, .375))
  expect_equal(aggregate_moments(original, order = 1:2)$value,
               aggregate_moments(paid, order = 1:2)$value)
  expect_equal(aggregate_moments(original, order = 2, type = "central")$value,
               aggregate_moments(paid, order = 2, type = "central")$value)
  expect_equal(original$count_basis, "loss")
})

test_that("registered binomial, negative-binomial and geometric thinning works", {
  severity <- apply_coverage(
    severity_spec("discrete", values = c(0, 10), probs = c(.2, .8)),
    coverage_spec())
  spec <- function(f) collective_risk_spec(f, severity)
  bin <- to_payment_basis(spec(frequency_spec("binomial", size = 10, prob = .5)))
  nb <- to_payment_basis(spec(frequency_spec("negative_binomial", size = 2, prob = .4)))
  geom <- to_payment_basis(spec(frequency_spec("geometric", prob = .4)))
  expect_equal(bin$frequency$params$prob, .4)
  expected <- .4 / (.4 + .6 * .8)
  expect_equal(nb$frequency$params$prob, expected)
  expect_equal(geom$frequency$params$prob, expected)
})

test_that("finite empirical frequency is thinned exactly", {
  freq <- empirical_frequency(c(0, 1, 2, 2))
  severity <- apply_coverage(
    severity_spec("discrete", values = c(0, 10), probs = c(.2, .8)),
    coverage_spec())
  paid <- to_payment_basis(collective_risk_spec(freq, severity))
  expect_identical(paid$frequency$dist, "empirical")
  expect_equal(distribution_moments(paid$frequency, 1)$value,
               .8 * distribution_moments(freq, 1)$value)
})

test_that("q=0 has no fictitious positive severity", {
  freq <- frequency_spec("poisson", lambda = 4)
  x <- severity_spec("discrete", values = c(25, 50), probs = c(.5, .5))
  model <- collective_risk_spec(freq,
            apply_coverage(x, coverage_spec(deductible = 100)))
  paid <- to_payment_basis(model)
  expect_s3_class(paid, "tidyact_zero_payment_basis_spec")
  expect_false(inherits(paid, "tidyact_collective_risk_spec"))
  expect_identical(paid$thinning$probability, 0)
  expect_null(paid$positive_severity)
  expect_identical(paid$metadata$aggregate_loss, 0)
})

test_that("q=1 preserves frequency and makes an executable collective model", {
  freq <- frequency_spec("poisson", lambda = 2)
  x <- severity_spec("gamma", shape = 2, rate = .01)
  model <- collective_risk_spec(freq, apply_coverage(x, coverage_spec()))
  paid <- to_payment_basis(model)
  expect_equal(paid$thinning$probability, 1)
  expect_equal(paid$frequency$params$lambda, 2)
  expect_equal(distribution_moments(paid$severity, 1)$value, 200)
  expect_equal(aggregate_moments(model, 1)$value,
               aggregate_moments(paid, 1)$value)
})

test_that("direct upper-tail probability survives cancellation", {
  model <- collective_risk_spec(frequency_spec("poisson", lambda = 4),
    apply_coverage(severity_spec("exponential", rate = 1),
                   coverage_spec(deductible = 40)))
  paid <- to_payment_basis(model)
  expect_equal(paid$thinning$probability, exp(-40), tolerance = 1e-13)
  expect_gt(paid$thinning$probability, 0)
})

test_that("basis checks and custom-frequency restrictions are explicit", {
  f <- frequency_spec("poisson", lambda = 1)
  severity <- apply_coverage(severity_spec("degenerate", value = 10),
                             coverage_spec())
  model <- collective_risk_spec(f, severity, count_basis = "payment")
  expect_error(to_payment_basis(model), "count_basis")
  custom <- frequency_spec("custom", custom = list(
    pmf = function(x, params, log = FALSE) stats::dpois(x, 1, log = log)
  ))
  model2 <- collective_risk_spec(custom, severity)
  # With q = 1 there is no actual thinning: the custom frequency is unchanged.
  expect_identical(to_payment_basis(model2)$frequency, custom)

  # With 0 < q < 1, custom frequency has no registered thinning adapter.
  severity_with_zeros <- apply_coverage(
    severity_spec("discrete", values = c(0, 10), probs = c(.2, .8)),
    coverage_spec()
  )
  model3 <- collective_risk_spec(custom, severity_with_zeros)
  expect_error(to_payment_basis(model3), "not yet registered")
})

test_that("conditioning removes the zero atom from uncapped continuous payments", {
  x <- severity_spec("exponential", rate = .01)
  model <- collective_risk_spec(frequency_spec("poisson", lambda = 5),
      apply_coverage(x, coverage_spec(deductible = 100)))
  paid <- to_payment_basis(model)
  expect_identical(paid$severity$support$type, "continuous")
  expect_equal(paid$thinning$probability, exp(-1))
  expect_equal(distribution_moments(paid$severity, 1)$value, 100)
  expect_equal(aggregate_moments(paid, order = 1)$value,
               aggregate_moments(model, order = 1)$value)
})

test_that("franchise cap below deductible gives a constant conditional payment", {
  x <- severity_spec("exponential", rate = .01)
  model <- collective_risk_spec(frequency_spec("poisson", lambda = 5),
    apply_coverage(x, coverage_spec(deductible = 100, limit = 50,
      deductible_type = "franchise", coinsurance = .8)))
  paid <- to_payment_basis(model)
  expect_identical(paid$severity$dist, "degenerate")
  expect_equal(paid$severity$params$value, 40)
  expect_equal(paid$thinning$probability, exp(-1))
})

test_that("finite-binomial aggregate PMFs coincide across payment bases", {
  freq <- frequency_spec("binomial", size = 3, prob = .5)
  loss <- severity_spec("discrete", values = c(0, 10),
                        probs = c(.2, .8))
  per_loss <- collective_risk_spec(freq,
                     apply_coverage(loss, coverage_spec(deductible = 5)))
  per_payment <- to_payment_basis(per_loss)
  a <- aggregate_distribution(per_loss, method = "convolution")
  b <- aggregate_distribution(per_payment, method = "convolution")
  expect_equal(a$grid$loss, b$grid$loss)
  expect_equal(a$grid$probability, b$grid$probability, tolerance = 1e-12)
  expect_equal(a$tail_probability, 0)
  expect_equal(b$tail_probability, 0)
})
test_that("thinning never silently rounds a positive payment count to zero", {
  paid <- apply_coverage(severity_spec("discrete", values = c(0, 1e20),
    probs = c(1, 1e-20)), coverage_spec())
  for (family in c("negative_binomial", "geometric")) {
    f <- if (family == "negative_binomial")
      frequency_spec(family, size = 2, prob = .5) else
      frequency_spec(family, prob = .5)
    model <- collective_risk_spec(f, paid)
    expect_equal(aggregate_moments(model, order = 1)$value,
                 if (family == "negative_binomial") 2 else 1)
    expect_error(to_payment_basis(model), "cannot be represented")
  }
})

test_that("thinning preserves representable tiny negative-binomial probabilities", {
  paid <- apply_coverage(severity_spec("discrete", values = c(0, 1),
    probs = c(.5, .5)), coverage_spec())
  for (family in c("negative_binomial", "geometric")) {
    f <- if (family == "negative_binomial")
      frequency_spec(family, size = 2, prob = 1e-310) else
      frequency_spec(family, prob = 1e-310)
    out <- to_payment_basis(collective_risk_spec(f, paid))
    expect_equal(out$frequency$params$prob / 1e-310, 2, tolerance = 1e-12)
  }
})
