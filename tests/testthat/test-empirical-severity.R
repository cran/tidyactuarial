test_that("empirical_severity builds exact empirical probabilities", {
  x <- empirical_severity(c(0, 500, 500, 1250.5, 3000))

  expect_s3_class(x, "tidyact_severity_spec")
  expect_identical(x$dist, "empirical")
  expect_equal(x$params$values, c(0, 500, 1250.5, 3000))
  expect_equal(x$params$probs, c(.2, .4, .2, .2))
  expect_equal(x$data$amount, x$params$values)
  expect_equal(x$data$probability, x$params$probs)
  expect_identical(x$support$type, "discrete")
})

test_that("weighted severity tables equal expanded data", {
  expanded <- empirical_severity(c(0, 500, 500, 1250.5, 3000))
  tab <- tibble::tibble(
    payment = c(0, 500, 1250.5, 3000),
    claims = c(1, 2, 1, 1)
  )
  weighted <- empirical_severity(tab, amount = payment, weight = claims)
  expect_equal(weighted$params, expanded$params)
})

test_that("zero weights are omitted without modifying positive observations", {
  x <- empirical_severity(c(0, 10, 100), weight = c(1, 0, 3))
  expect_equal(x$params$values, c(0, 100))
  expect_equal(x$params$probs, c(.25, .75))
})

test_that("empirical_severity validates amounts", {
  expect_error(empirical_severity(numeric()), "nonempty")
  expect_error(empirical_severity(c(0, -1)), "nonnegative")
  expect_error(empirical_severity(c(0, NA_real_)), "nonnegative")
  expect_error(empirical_severity(c(0, Inf)), "finite")
  expect_error(empirical_severity(c("a", "b")), "numeric")
})

test_that("empirical_severity validates weights", {
  expect_error(empirical_severity(c(10, 20), weight = 1), "one finite")
  expect_error(empirical_severity(c(10, 20), weight = c(1, -1)), "nonnegative")
  expect_error(empirical_severity(c(10, 20), weight = c(0, 0)), "positive weight")
  expect_error(empirical_severity(c(10, 20), weight = c(1, Inf)), "finite")
})

test_that("data-frame severity input requires an amount column", {
  d <- tibble::tibble(x = c(10, 20))
  expect_error(empirical_severity(d), "amount")
  expect_error(empirical_severity(c(10, 20), amount = x), "omitted")
})

test_that("empirical severity PMF, CDF and quantiles are coherent", {
  x <- empirical_severity(c(0, 100, 100, 300))
  entry <- tidyactuarial:::.severity_registry()$empirical

  expect_equal(entry$pmf(c(0, 50, 100, 300), x$params), c(.25, 0, .5, .25))
  expect_equal(entry$cdf(c(-1, 0, 99, 100, 300), x$params), c(0, .25, .25, .75, 1))
  expect_equal(entry$quantile(c(0, .25, .26, .75, .76, 1), x$params), c(0, 0, 100, 100, 300, 300))
  expect_equal(
    entry$quantile(log(c(.25, .75)), x$params, log.p = TRUE),
    c(0, 100)
  )
})

test_that("empirical severity analytical operations match direct sums", {
  x <- empirical_severity(c(0, 100, 100, 300))
  entry <- tidyactuarial:::.severity_registry()$empirical

  expect_equal(entry$raw_moment(1, x$params), 125)
  expect_equal(entry$raw_moment(2, x$params), 27500)
  expect_equal(
    entry$limited_moment(1, 150, x$params),
    mean(pmin(c(0, 100, 100, 300), 150))
  )
  expect_equal(entry$mgf(0, x$params), 1)
  expect_equal(entry$laplace(0, x$params), 1)
})

test_that("empirical severity RNG samples only observed support", {
  x <- empirical_severity(c(0, 100, 100, 300))
  entry <- tidyactuarial:::.severity_registry()$empirical
  set.seed(444)
  draws <- entry$rng(100, x$params)
  expect_length(draws, 100)
  expect_true(all(draws %in% x$params$values))
})

test_that("empirical severity keeps decimal amounts distinct", {
  x <- empirical_severity(c(10, 10.0000001, 10))
  expect_length(x$params$values, 2)
  expect_equal(x$params$probs, c(2 / 3, 1 / 3))
})
