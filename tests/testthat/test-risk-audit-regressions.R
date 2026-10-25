test_that("aggregate engines preserve monetary scale and truncation", {
  for (scale in c(1, 1e-16, 1e16)) {
    model <- collective_risk_spec(
      frequency_spec("binomial", size = 2, prob = 0.5),
      severity_spec("discrete", values = c(1, 2) * scale,
                    probs = c(0.5, 0.5)))
    for (method in c("convolution", "panjer", "fft")) {
      out <- aggregate_distribution(model, method = method)
      # Direct expansion of (.5 + .25 z + .25 z^2)^2.
      expect_equal(out$grid$loss / scale, 0:4)
      expect_equal(out$grid$probability, c(.25, .25, .3125, .125, .0625))
      cut <- aggregate_distribution(model, method = method,
                                  max_loss = 2 * scale)
      expect_equal(cut$grid$loss / scale, 0:2)
      expect_equal(cut$tail_probability, .1875)
    }
  }
})

test_that("rare positive masses retain their expected-loss contribution", {
  model <- collective_risk_spec(
    frequency_spec("binomial", size = 1, prob = 1),
    severity_spec("discrete", values = c(0, 1e20), probs = c(1, 1e-20)))
  out <- aggregate_distribution(model, method = "convolution")
  expect_equal(sum(out$grid$loss * out$grid$probability), 1)
  expect_equal(stop_loss_premium(out, 5e19)$premium, .5)
})

test_that("exponential coverage stop-loss respects high limits and no limit", {
  for (limit in c(.5, 2, 1e15, Inf)) {
    for (kind in c("ordinary", "franchise")) {
      paid <- apply_coverage(severity_spec("exponential", rate = 1),
        coverage_spec(deductible = 1, limit = limit,
                      coinsurance = .5, deductible_type = kind))
      retentions <- c(0, .25, .75, 2)
      expected <- vapply(retentions, function(d) {
        # Independent integration on the ground-up exponential variable.
        stats::integrate(function(z) {
          y <- if (kind == "ordinary") .5 * pmin(pmax(z - 1, 0), limit)
            else .5 * pmin(z, limit) * (z > 1)
          pmax(y - d, 0) * exp(-z)
        }, 0, Inf, subdivisions = 1000, rel.tol = 1e-9)$value
      }, numeric(1))
      expect_equal(stop_loss_premium(paid, retentions)$premium, expected,
                   tolerance = 1e-8)
    }
  }
})

test_that("covered exponential stop-loss scales with inflation and coinsurance", {
  paid <- apply_coverage(severity_spec("exponential", rate = .7),
    coverage_spec(deductible = 2, limit = 3, inflation_factor = 1.5,
                  coinsurance = .8, index_deductible = TRUE,
                  index_limit = TRUE))
  # Effective Z has rate .7/1.5, deductible 3, layer limit 4.5.
  r <- .7 / 1.5
  expected <- .8 / r * exp(-r * (3 + 1 / .8)) *
    (-expm1(-r * (4.5 - 1 / .8)))
  expect_equal(stop_loss_premium(paid, 1)$premium, expected)
  expect_equal(tail_value_at_risk(paid, .99)$tail_value_at_risk, 3.6)
})
