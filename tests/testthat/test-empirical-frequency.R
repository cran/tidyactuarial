test_that("empirical_frequency builds exact empirical probabilities", {
  x <- empirical_frequency(c(0, 0, 0, 1, 1, 2))

  expect_s3_class(x, "tidyact_frequency_spec")
  expect_identical(x$dist, "empirical")
  expect_equal(x$params$values, c(0, 1, 2))
  expect_equal(x$params$probs, c(.5, 1 / 3, 1 / 6))
  expect_equal(x$data$count, c(0, 1, 2))
  expect_equal(x$data$weight, c(3, 2, 1))
  expect_equal(x$data$probability, c(.5, 1 / 3, 1 / 6))
  expect_identical(x$support$type, "count")
})

test_that("weighted tables and expanded observations give the same frequency", {
  expanded <- empirical_frequency(c(0, 0, 0, 1, 1, 2))
  tab <- tibble::tibble(claims = c(0, 1, 2), policies = c(3, 2, 1))
  weighted <- empirical_frequency(tab, count = claims, weight = policies)

  expect_equal(weighted$params, expanded$params)
})

test_that("zero-weight count rows are omitted", {
  x <- empirical_frequency(c(0, 1, 5), weight = c(2, 1, 0))
  expect_equal(x$params$values, c(0, 1))
  expect_equal(x$params$probs, c(2 / 3, 1 / 3))
  expect_equal(x$support$upper, 1)
})

test_that("empirical_frequency validates counts", {
  expect_error(empirical_frequency(numeric()), "nonempty")
  expect_error(empirical_frequency(c(0, -1)), "nonnegative integers")
  expect_error(empirical_frequency(c(0, 1.5)), "nonnegative integers")
  expect_error(empirical_frequency(c(0, NA_real_)), "nonnegative integers")
  expect_error(empirical_frequency(c(0, Inf)), "nonnegative integers")
  expect_error(empirical_frequency(c(0, 2^53 + 2)), "nonnegative integers")
})

test_that("empirical_frequency validates weights", {
  expect_error(empirical_frequency(c(0, 1), weight = 1), "one finite")
  expect_error(empirical_frequency(c(0, 1), weight = c(1, -1)), "nonnegative")
  expect_error(empirical_frequency(c(0, 1), weight = c(0, 0)), "positive weight")
  expect_error(empirical_frequency(c(0, 1), weight = c(1, Inf)), "finite")
})

test_that("data-frame frequency input requires a count column", {
  d <- tibble::tibble(n = c(0, 1, 2))
  expect_error(empirical_frequency(d), "count")
  expect_error(empirical_frequency(c(0, 1), count = n), "omitted")
})

test_that("empirical frequency PMF, CDF and quantiles are coherent", {
  x <- empirical_frequency(c(0, 0, 1, 3))
  entry <- tidyactuarial:::.frequency_registry()$empirical

  expect_true(is.function(entry$pmf))
  expect_equal(entry$pmf(c(0, 1, 2, 3), x$params), c(.5, .25, 0, .25))
  expect_equal(entry$cdf(c(-1, 0, 1, 2, 3), x$params), c(0, .5, .75, .75, 1))
  expect_equal(entry$quantile(c(0, .5, .51, .75, .76, 1), x$params), c(0, 0, 1, 1, 3, 3))
  expect_equal(
    entry$quantile(log(c(.5, .75)), x$params, log.p = TRUE),
    c(0, 1)
  )
})

test_that("empirical frequency moments and PGF agree with direct definitions", {
  x <- empirical_frequency(c(0, 0, 1, 3))
  entry <- tidyactuarial:::.frequency_registry()$empirical

  expect_equal(entry$raw_moment(1, x$params), 1)
  expect_equal(entry$raw_moment(2, x$params), 2.5)
  expect_equal(entry$factorial_moment(1, x$params), 1)
  expect_equal(entry$factorial_moment(2, x$params), 1.5)
  expect_equal(entry$pgf(1, x$params), 1)

  z <- .7
  expect_equal(
    entry$pgf(z, x$params),
    sum(x$params$probs * z^x$params$values)
  )
  expect_null(entry$panjer)
})

test_that("empirical frequency RNG samples only observed support", {
  x <- empirical_frequency(c(0, 0, 2, 5))
  entry <- tidyactuarial:::.frequency_registry()$empirical
  set.seed(321)
  draws <- entry$rng(100, x$params)
  expect_length(draws, 100)
  expect_true(all(draws %in% x$params$values))
})
