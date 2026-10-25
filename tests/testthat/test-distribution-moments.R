test_that("distribution_moments returns the stable tibble contract", {
  x <- frequency_spec("poisson", lambda = 4)
  out <- distribution_moments(x, order = c(3, 1, 3), type = "raw")

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("type", "order", "value", "status"))
  expect_identical(out$type, rep("raw", 3))
  expect_identical(out$order, c(3L, 1L, 3L))
  expect_identical(out$status, rep("finite", 3))
})

test_that("Poisson raw, factorial, central moments and cumulants are correct", {
  x <- frequency_spec("poisson", lambda = 4)

  raw <- distribution_moments(x, 0:4, "raw")
  expect_equal(raw$value, c(1, 4, 20, 116, 756))

  factorial <- distribution_moments(x, 0:4, "factorial")
  expect_equal(factorial$value, 4^(0:4))

  central <- distribution_moments(x, 0:4, "central")
  expect_equal(central$value, c(1, 0, 4, 4, 4 + 3 * 4^2), tolerance = 1e-12)

  cumulant <- distribution_moments(x, 1:4, "cumulant")
  expect_equal(cumulant$value, rep(4, 4), tolerance = 1e-12)
})

test_that("binomial moments recover mean and variance", {
  x <- frequency_spec("binomial", size = 20, prob = .3)
  raw <- distribution_moments(x, 1:2)
  mean <- 20 * .3
  variance <- 20 * .3 * .7

  expect_equal(raw$value[1], mean, tolerance = 1e-12)
  expect_equal(raw$value[2] - raw$value[1]^2, variance, tolerance = 1e-12)
})

test_that("severity raw moments agree with analytical identities", {
  expo <- severity_spec("exponential", rate = .2)
  expect_equal(
    distribution_moments(expo, 0:4)$value,
    factorial(0:4) / .2^(0:4),
    tolerance = 1e-10
  )

  gamma <- severity_spec("gamma", shape = 3, rate = .5)
  expected <- vapply(0:4, function(k) {
    exp(lgamma(3 + k) - lgamma(3) - k * log(.5))
  }, numeric(1))
  expect_equal(distribution_moments(gamma, 0:4)$value, expected, tolerance = 1e-10)
})

test_that("central moments of a deterministic severity are zero", {
  x <- severity_spec("degenerate", value = 2500)
  out <- distribution_moments(x, 0:5, type = "central")
  expect_equal(out$value, c(1, rep(0, 5)))
  expect_identical(out$status, rep("finite", 6))
})

test_that("finite discrete central moments are calculated directly", {
  x <- severity_spec("discrete", values = c(0, 10, 20), probs = c(.2, .5, .3))
  mu <- sum(x$params$values * x$params$probs)
  expected <- vapply(0:4, function(k) {
    if (k == 0) 1 else sum(x$params$probs * (x$params$values - mu)^k)
  }, numeric(1))
  expect_equal(distribution_moments(x, 0:4, "central")$value, expected)
})

test_that("heavy-tail moment existence is reported rather than approximated", {
  x <- severity_spec("pareto", shape = 1.5, scale = 2000)
  raw <- distribution_moments(x, 1:4)

  expect_identical(raw$status, c("finite", "infinite", "infinite", "infinite"))
  expect_true(is.finite(raw$value[1]))
  expect_true(all(is.infinite(raw$value[2:4])))

  central <- distribution_moments(x, 2:4, type = "central")
  expect_identical(central$status, rep("infinite", 3))
})

test_that("central moments are undefined when the mean does not exist", {
  x <- severity_spec("pareto", shape = .8, scale = 1000)
  out <- distribution_moments(x, 1:3, type = "central")
  expect_identical(out$status, rep("undefined", 3))
  expect_true(all(is.na(out$value)))
})

test_that("cumulants requiring divergent raw moments are undefined", {
  x <- severity_spec("pareto", shape = 1.5, scale = 1000)
  out <- distribution_moments(x, 1:3, type = "cumulant")
  expect_identical(out$status[1], "finite")
  expect_identical(out$status[2:3], c("undefined", "undefined"))
})

test_that("empirical distributions integrate with distribution_moments", {
  n <- empirical_frequency(c(0, 0, 1, 2))
  raw_n <- distribution_moments(n, 1:2)
  expect_equal(raw_n$value, c(.75, 1.25))

  fac_n <- distribution_moments(n, 1:2, "factorial")
  expect_equal(fac_n$value, c(.75, .5))

  x <- empirical_severity(c(0, 100, 100, 300))
  raw_x <- distribution_moments(x, 1:2)
  expect_equal(raw_x$value, c(125, 27500))

  central_x <- distribution_moments(x, 2, "central")
  expect_equal(central_x$value, mean((c(0, 100, 100, 300) - 125)^2))
})

test_that("factorial moments are restricted to frequency specifications", {
  x <- severity_spec("gamma", shape = 2, rate = 1)
  expect_error(distribution_moments(x, 1:2, "factorial"), "frequency")
})

test_that("distribution_moments validates x, order and type", {
  expect_error(distribution_moments(1), "frequency or severity")

  x <- frequency_spec("poisson", lambda = 2)
  expect_error(distribution_moments(x, numeric()), "nonempty")
  expect_error(distribution_moments(x, -1), "nonnegative")
  expect_error(distribution_moments(x, 1.5), "nonnegative")
  expect_error(distribution_moments(x, NA_real_), "nonnegative")
  expect_error(distribution_moments(x, 1, "variance"), "type")
  expect_error(distribution_moments(x, 0, "cumulant"), "positive")
})

test_that("custom distributions do not trigger hidden numerical integration", {
  x <- severity_spec(
    "custom",
    custom = list(
      pdf = function(x, params, log = FALSE) stats::dexp(x, 1, log = log),
      support = list(type = "continuous", lower = 0, upper = Inf)
    )
  )

  out <- distribution_moments(x, 1:3)
  expect_identical(out$status, rep("undefined", 3))
  expect_true(all(is.na(out$value)))
})

test_that("custom analytical moments are used when explicitly supplied", {
  x <- severity_spec(
    "custom",
    rate = 2,
    custom = list(
      pdf = function(x, params, log = FALSE) stats::dexp(x, params$rate, log = log),
      raw_moment = function(order, params) factorial(order) / params$rate^order,
      support = list(type = "continuous", lower = 0, upper = Inf)
    )
  )

  expect_equal(
    distribution_moments(x, 1:4)$value,
    factorial(1:4) / 2^(1:4)
  )
})
