test_that("severity_spec constructs built-in continuous severities", {
  cases <- list(
    exponential = list(rate = 0.002),
    gamma = list(shape = 2.5, rate = 0.001),
    lognormal = list(meanlog = 8, sdlog = 0.8),
    weibull = list(shape = 1.7, scale = 5000),
    pareto = list(shape = 3, scale = 2000),
    burr = list(shape1 = 1.5, shape2 = 2.5, scale = 1500),
    loglogistic = list(shape = 2.5, scale = 1500)
  )

  for (dist in names(cases)) {
    x <- do.call(severity_spec, c(list(dist = dist), cases[[dist]]))
    expect_s3_class(x, "tidyact_severity_spec")
    expect_s3_class(x, "tidyact_distribution_spec")
    expect_identical(x$schema_version, 1L)
    expect_identical(x$dist, dist)
    expect_identical(x$support$type, "continuous")
    expect_equal(x$support$lower, 0)
    expect_equal(x$support$upper, Inf)
    expect_null(x$data)
    expect_null(x$custom)
  }
})

test_that("severity_spec constructs degenerate and finite discrete severities", {
  fixed <- severity_spec("degenerate", value = 1000)
  expect_identical(fixed$support$type, "discrete")
  expect_equal(fixed$support$lower, 1000)
  expect_equal(fixed$support$upper, 1000)

  d <- severity_spec(
    "discrete",
    values = c(0, 500, 2500),
    probs = c(0.1, 0.6, 0.3)
  )
  expect_identical(d$support$type, "discrete")
  expect_equal(d$support$lower, 0)
  expect_equal(d$support$upper, 2500)
})

test_that("severity distribution names are normalized and labels are stored", {
  x <- severity_spec("  GaMmA  ", shape = 2, rate = 0.1, label = "Losses")
  expect_identical(x$dist, "gamma")
  expect_identical(x$label, "Losses")
})

test_that("severity_spec rejects malformed distribution identifiers", {
  expect_error(severity_spec(NULL), "dist")
  expect_error(severity_spec(""), "dist")
  expect_error(severity_spec(NA_character_), "dist")
  expect_error(severity_spec(c("gamma", "weibull")), "dist")
  expect_error(severity_spec("normal", mean = 0, sd = 1), "must be one of")
})

test_that("built-in severities reject missing, extra and invalid parameters", {
  expect_error(severity_spec("gamma", shape = 2), "Missing")
  expect_error(severity_spec("gamma", shape = 2, rate = 1, scale = 1), "Unused")
  expect_error(severity_spec("exponential", rate = 0), "rate")
  expect_error(severity_spec("gamma", shape = -1, rate = 1), "shape")
  expect_error(severity_spec("gamma", shape = 1, rate = Inf), "rate")
  expect_error(severity_spec("lognormal", meanlog = 0, sdlog = 0), "sdlog")
  expect_error(severity_spec("weibull", shape = 0, scale = 1), "shape")
  expect_error(severity_spec("pareto", shape = 1, scale = 0), "scale")
  expect_error(severity_spec("burr", shape1 = 1, shape2 = 0, scale = 1), "shape2")
  expect_error(severity_spec("loglogistic", shape = 1, scale = -1), "scale")
  expect_error(severity_spec("degenerate", value = -1), "value")
})

test_that("finite discrete severity validation is strict", {
  expect_error(
    severity_spec("discrete", values = c(0, 2, 1), probs = c(.2, .3, .5)),
    "strictly increasing"
  )
  expect_error(
    severity_spec("discrete", values = c(0, 1), probs = c(.5, .4)),
    "sum to one"
  )
  expect_error(
    severity_spec("discrete", values = c(0, 1), probs = c(1, 0)),
    "strictly positive"
  )
  expect_error(
    severity_spec("discrete", values = c(0, Inf), probs = c(.5, .5)),
    "finite nonnegative"
  )
})

test_that("stats-backed severity operations agree with stats", {
  registry <- tidyactuarial:::.severity_registry()
  x <- c(0, 100, 500, 1000)
  p <- c(.01, .25, .5, .9)

  pars <- list(rate = .002)
  expect_equal(registry$exponential$pdf(x, pars), stats::dexp(x, rate = .002))
  expect_equal(registry$exponential$cdf(x, pars), stats::pexp(x, rate = .002))
  expect_equal(registry$exponential$quantile(p, pars), stats::qexp(p, rate = .002))

  pars <- list(shape = 2.5, rate = .001)
  expect_equal(registry$gamma$pdf(x, pars), stats::dgamma(x, 2.5, rate = .001))
  expect_equal(registry$gamma$cdf(x, pars), stats::pgamma(x, 2.5, rate = .001))
  expect_equal(registry$gamma$quantile(p, pars), stats::qgamma(p, 2.5, rate = .001))

  pars <- list(meanlog = 7, sdlog = .6)
  expect_equal(registry$lognormal$pdf(x, pars), stats::dlnorm(x, 7, .6))
  expect_equal(registry$lognormal$cdf(x, pars), stats::plnorm(x, 7, .6))
  expect_equal(registry$lognormal$quantile(p, pars), stats::qlnorm(p, 7, .6))

  pars <- list(shape = 1.8, scale = 1000)
  expect_equal(registry$weibull$pdf(x, pars), stats::dweibull(x, 1.8, 1000))
  expect_equal(registry$weibull$cdf(x, pars), stats::pweibull(x, 1.8, 1000))
  expect_equal(registry$weibull$quantile(p, pars), stats::qweibull(p, 1.8, 1000))
})

test_that("Pareto II, Burr XII and log-logistic quantiles invert their CDFs", {
  registry <- tidyactuarial:::.severity_registry()
  probs <- c(.01, .2, .5, .9, .99)

  cases <- list(
    pareto = list(shape = 3.2, scale = 1200),
    burr = list(shape1 = 1.7, shape2 = 2.3, scale = 1500),
    loglogistic = list(shape = 2.5, scale = 1800)
  )

  for (dist in names(cases)) {
    q <- registry[[dist]]$quantile(probs, cases[[dist]])
    expect_equal(
      registry[[dist]]$cdf(q, cases[[dist]]),
      probs,
      tolerance = 1e-10
    )
  }
})

test_that("discrete severity probability functions are coherent", {
  entry <- tidyactuarial:::.severity_registry()$discrete
  pars <- list(values = c(0, 10, 25), probs = c(.2, .5, .3))

  expect_equal(entry$pmf(c(0, 5, 10, 25), pars), c(.2, 0, .5, .3))
  expect_equal(entry$cdf(c(-1, 0, 9, 10, 25), pars), c(0, .2, .2, .7, 1))
  expect_equal(entry$quantile(c(0, .2, .21, .7, 1), pars), c(0, 0, 10, 10, 25))
  expect_equal(
    entry$quantile(log(c(.2, .7)), pars, log.p = TRUE),
    c(0, 10)
  )
})

test_that("degenerate severity is mathematically degenerate", {
  entry <- tidyactuarial:::.severity_registry()$degenerate
  pars <- list(value = 100)

  expect_equal(entry$pmf(c(99, 100, 101), pars), c(0, 1, 0))
  expect_equal(entry$cdf(c(99, 100, 101), pars), c(0, 1, 1))
  expect_equal(entry$raw_moment(3, pars), 100^3)
  expect_equal(entry$limited_moment(2, 40, pars), 40^2)
  expect_equal(entry$laplace(0, pars), 1)
})

test_that("custom severity stores only explicitly supplied operations", {
  custom <- list(
    pdf = function(x, params, log = FALSE) stats::dexp(x, params$rate, log = log),
    cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
      stats::pexp(q, params$rate, lower.tail = lower.tail, log.p = log.p)
    },
    raw_moment = function(order, params) factorial(order) / params$rate^order,
    support = list(type = "continuous", lower = 0, upper = Inf)
  )

  x <- severity_spec("custom", rate = .1, custom = custom)
  expect_identical(x$dist, "custom")
  expect_true(is.function(x$custom$pdf))
  expect_true(is.function(x$custom$raw_moment))
  expect_null(x$custom$support)
  expect_identical(x$support$type, "continuous")
})

test_that("custom definitions are validated", {
  expect_error(severity_spec("custom", custom = NULL), "custom")
  expect_error(
    severity_spec("custom", custom = list(pdf = function(...) 1)),
    "support"
  )
  expect_error(
    severity_spec(
      "custom",
      custom = list(pdf = 1, support = list(type = "continuous", lower = 0, upper = Inf))
    ),
    "function"
  )
  expect_error(
    severity_spec(
      "gamma", shape = 2, rate = 1,
      custom = list(pdf = function(...) 1)
    ),
    "custom"
  )
})

test_that("severity print method is informative and returns invisibly", {
  x <- severity_spec("gamma", shape = 2, rate = .1, label = "Claims")
  out <- capture.output(returned <- print(x))
  expect_identical(returned, x)
  expect_true(any(grepl("tidyact_severity_spec", out, fixed = TRUE)))
  expect_true(any(grepl("Distribution: Gamma", out, fixed = TRUE)))
  expect_true(any(grepl("Label: Claims", out, fixed = TRUE)))
})
