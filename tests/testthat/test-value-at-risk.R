test_that("VaR uses a left quantile with atoms, not linear interpolation", {
  loss <- c(0, 100, 100, 300)
  out <- value_at_risk(loss, level = c(0.5, 0.75, 0.9))
  expect_equal(out$value_at_risk, c(100, 100, 300))
  expect_equal(out$level, c(0.5, 0.75, 0.9))
})

test_that("VaR accepts simulations and analytic continuous severities", {
  sims <- simulate_aggregate(collective_risk_spec(
    frequency_spec("binomial", size = 1, prob = 1),
    severity_spec("degenerate", value = 50)
  ), 30, seed = 8)
  expect_equal(value_at_risk(sims, level = 0.95)$value_at_risk, 50)
  x <- severity_spec("exponential", rate = 0.01)
  expect_equal(value_at_risk(x, 0.95)$value_at_risk,
               stats::qexp(0.95, rate = 0.01))
})

test_that("VaR does not renormalize unresolved aggregate tails", {
  grid <- structure(list(grid = tibble::tibble(
    loss = c(0, 100), probability = c(0.5, 0.45)),
    tail_probability = 0.05), class = "tidyact_aggregate_distribution")
  expect_equal(value_at_risk(grid, 0.90)$value_at_risk, 100)
  expect_error(value_at_risk(grid, 0.98), "uncomputed aggregate tail")
})

test_that("VaR refuses invalid levels and missing simulation columns", {
  expect_error(value_at_risk(c(0, 1), level = 1), "level")
  expect_error(value_at_risk(c(0, 1), level = -0.1), "level")
  expect_error(value_at_risk(tibble::tibble(a = 1:3)), "value")
})

test_that("VaR respects strict probability boundaries without an epsilon shift", {
  x <- severity_spec("discrete", values = c(0, 100), probs = c(.5, .5))
  expect_equal(value_at_risk(x, c(.5, .5 + .Machine$double.eps))$value_at_risk,
               c(0, 100))
  rare <- severity_spec("discrete", values = c(1, 2), probs = c(1e-20, 1))
  expect_equal(value_at_risk(rare, c(1e-20, 2e-20))$value_at_risk, c(1, 2))
})

test_that("VaR rejects levels immediately beyond retained aggregate mass", {
  grid <- structure(list(grid = tibble::tibble(
    loss = c(0, 100), probability = c(.5, .45)),
    tail_probability = .05), class = "tidyact_aggregate_distribution")
  retained <- sum(grid$grid$probability)
  expect_equal(value_at_risk(grid, retained)$value_at_risk, 100)
  expect_error(value_at_risk(grid, retained + .Machine$double.eps),
               "uncomputed aggregate tail")
})

test_that("FFT count truncation cannot silently move a loss quantile", {
  model <- collective_risk_spec(frequency_spec("poisson", lambda = 2),
    severity_spec("discrete", values = c(0, 1), probs = c(.9, .1)))
  coarse <- aggregate_distribution(model, method = "fft", tol = .1)
  # Bernoulli thinning gives S ~ Poisson(.2) exactly. Omitted large
  # counts still contribute at zero; they are not exclusively large losses.
  expect_equal(stats::qpois(.8, .2), 0)
  expect_error(value_at_risk(coarse, .8), "count truncation")
  # Where the retained mass and omitted-mass bound select the same atom,
  # the quantile remains identifiable despite count truncation.
  expect_equal(value_at_risk(coarse, c(.7, .9))$value_at_risk,
               stats::qpois(c(.7, .9), .2))
  fine <- aggregate_distribution(model, method = "fft", tol = 1e-9)
  expect_equal(value_at_risk(fine, .8)$value_at_risk, stats::qpois(.8, .2))
})

test_that("VaR uses the minimum possible loss of omitted FFT counts", {
  model <- collective_risk_spec(frequency_spec("poisson", lambda = 2),
    severity_spec("discrete", values = c(1, 2), probs = c(.5, .5)))
  out <- aggregate_distribution(model, method = "fft", tol = .1)
  lower <- out$metadata$count_truncation + 1
  below <- sum(out$grid$probability[out$grid$loss < lower])
  at <- sum(out$grid$probability[out$grid$loss == lower])
  p <- below + min(at, out$metadata$count_tail_probability) / 2
  # All omitted paths have at least 'lower' unit claims. They can add
  # probability at this amount, but cannot move its quantile below it.
  expect_equal(value_at_risk(out, p)$value_at_risk, lower)
})
