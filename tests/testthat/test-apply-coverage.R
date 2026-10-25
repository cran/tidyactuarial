test_that("finite discrete coverage uses the installed discrete provider exactly", {
  x <- severity_spec("discrete", values = c(50, 150, 250),
                     probs = c(0.2, 0.5, 0.3))
  cov <- coverage_spec(deductible = 100)
  result <- apply_coverage(x, cov)
  expect_s3_class(result, "tidyact_payment_spec")
  expect_identical(result$basis, "per_loss")
  expect_identical(result$dist, "discrete")
  expect_identical(result$support$type, "discrete")
  expect_equal(result$params$values, c(0, 50, 150))
  expect_equal(result$params$probs, c(0.2, 0.5, 0.3))
  expect_equal(distribution_moments(result, order = 1)$value, 70)
  expect_equal(distribution_moments(result, order = 2)$value, 8000)
  expect_s3_class(collective_risk_spec(frequency_spec("poisson", lambda = 3),
                                        result), "tidyact_collective_risk_spec")
})

test_that("discrete franchise retains threshold mass and cap mass exactly", {
  x <- severity_spec("discrete", values = c(100, 150, 200, 350),
                     probs = c(0.1, 0.2, 0.3, 0.4))
  y <- apply_coverage(x, coverage_spec(deductible = 150, limit = 200,
                        deductible_type = "franchise", coinsurance = 0.5))
  expect_equal(y$params$values, c(0, 100))
  expect_equal(y$params$probs, c(0.3, 0.7))
})

test_that("continuous exponential payment is mixed with exact CDF and atoms", {
  x <- severity_spec("exponential", rate = 0.01)
  y <- apply_coverage(x, coverage_spec(deductible = 100, limit = 150))
  expect_identical(y$dist, "custom")
  expect_identical(y$support$type, "mixed")
  cdf <- y$custom$cdf
  pmf <- y$custom$pmf
  expect_equal(cdf(c(-1, 0, 50, 149, 150, Inf), y$params),
               c(0, stats::pexp(100, rate = .01),
                 stats::pexp(150, rate = .01),
                 stats::pexp(249, rate = .01), 1, 1))
  expect_equal(pmf(c(0, 50, 150), y$params),
               c(stats::pexp(100, rate = .01), 0,
                 stats::pexp(250, rate = .01, lower.tail = FALSE)))
  expect_equal(distribution_moments(y, order = 1)$value,
               exp(-1) / 0.01 * stats::pgamma(150, shape = 1, rate = .01))
})

test_that("franchise has exact zero mass, upper atom and positive gap", {
  x <- severity_spec("exponential", rate = .01)
  y <- apply_coverage(x, coverage_spec(deductible = 100,
                        limit = 150, deductible_type = "franchise"))
  expect_equal(y$custom$cdf(50, y$params), stats::pexp(100, rate = .01))
  expect_equal(y$custom$pdf(c(50, 100, 125, 150), y$params),
               c(0, 0, stats::dexp(125, rate = .01), 0))
  expect_equal(y$custom$pmf(150, y$params), exp(-1.5))
})

test_that("zero coverage and full coverage use valid registered distributions", {
  x <- severity_spec("gamma", shape = 2, rate = .01)
  y0 <- apply_coverage(x, coverage_spec(limit = 0))
  expect_equal(y0$params$values, 0)
  expect_equal(y0$params$probs, 1)
  expect_equal(distribution_moments(y0, 1)$value, 0)
  y <- apply_coverage(x, coverage_spec())
  expect_identical(y$support$type, "continuous")
  expect_equal(y$custom$cdf(100, y$params), stats::pgamma(100, 2, rate = .01))
  expect_equal(distribution_moments(y, 1)$value, 200)
})

test_that("unsupported exact source provider fails without silent approximation", {
  bad <- severity_spec("custom", custom = list(
    support = list(type = "continuous", lower = 0, upper = Inf),
    pdf = function(x, params, log = FALSE) stats::dexp(x, log = log)
  ))
  expect_error(apply_coverage(bad, coverage_spec(deductible = 100)),
               "exact CDF, quantile and RNG")
})
