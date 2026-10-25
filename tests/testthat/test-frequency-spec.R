# Tests for frequency_spec -------------------------------------------------

# Helpers -----------------------------------------------------------------

frequency_registry <- function() {
  tidyactuarial:::.frequency_registry()
}


expect_valid_frequency_spec <- function(x, dist) {
  expect_s3_class(x, "tidyact_frequency_spec")
  expect_s3_class(x, "tidyact_distribution_spec")

  expect_identical(x$schema_version, 1L)
  expect_identical(x$dist, dist)

  expect_true(is.list(x$params))
  expect_true(is.list(x$support))
  expect_null(x$data)
  expect_true(is.list(x$metadata))

  expect_named(
    x,
    c(
      "schema_version",
      "dist",
      "params",
      "support",
      "data",
      "custom",
      "label",
      "metadata"
    ),
    ignore.order = FALSE
  )
}


expect_count_support <- function(
    x,
    lower = 0,
    upper = Inf,
    step = 1
) {
  expect_identical(x$support$type, "count")
  expect_equal(x$support$lower, lower)
  expect_equal(x$support$upper, upper)
  expect_equal(x$support$step, step)
}


expect_panjer_recursion <- function(
    entry,
    params,
    k
) {
  panjer <- entry$panjer(params)

  expect_true(is.list(panjer))

  expect_named(
    panjer,
    c("a", "b", "p0", "p1"),
    ignore.order = FALSE
  )

  expect_equal(
    panjer$p0,
    entry$pmf(0, params)
  )

  expect_equal(
    panjer$p1,
    entry$pmf(1, params)
  )

  lhs <- entry$pmf(
    k,
    params
  )

  rhs <- (
    panjer$a +
      panjer$b / k
  ) *
    entry$pmf(
      k - 1,
      params
    )

  expect_equal(
    lhs,
    rhs,
    tolerance = 1e-12
  )
}


# Constructor --------------------------------------------------------------

test_that("frequency_spec constructs Poisson specifications", {
  x <- frequency_spec(
    "poisson",
    lambda = 4
  )

  expect_valid_frequency_spec(
    x,
    "poisson"
  )

  expect_equal(
    x$params,
    list(lambda = 4)
  )

  expect_count_support(
    x,
    lower = 0,
    upper = Inf,
    step = 1
  )

  expect_null(x$custom)
  expect_null(x$label)
})


test_that("frequency_spec constructs binomial specifications", {
  x <- frequency_spec(
    "binomial",
    size = 20,
    prob = 0.25
  )

  expect_valid_frequency_spec(
    x,
    "binomial"
  )

  expect_equal(
    x$params,
    list(
      size = 20,
      prob = 0.25
    )
  )

  expect_count_support(
    x,
    lower = 0,
    upper = 20,
    step = 1
  )
})


test_that("frequency_spec constructs negative binomial specifications", {
  x <- frequency_spec(
    "negative_binomial",
    size = 2.5,
    prob = 0.4
  )

  expect_valid_frequency_spec(
    x,
    "negative_binomial"
  )

  expect_equal(
    x$params,
    list(
      size = 2.5,
      prob = 0.4
    )
  )

  expect_count_support(
    x,
    lower = 0,
    upper = Inf,
    step = 1
  )
})


test_that("frequency_spec constructs geometric specifications", {
  x <- frequency_spec(
    "geometric",
    prob = 0.25
  )

  expect_valid_frequency_spec(
    x,
    "geometric"
  )

  expect_equal(
    x$params,
    list(prob = 0.25)
  )

  expect_count_support(
    x,
    lower = 0,
    upper = Inf,
    step = 1
  )
})


test_that("distribution identifiers are normalized", {
  x <- frequency_spec(
    "  PoIsSoN  ",
    lambda = 3
  )

  expect_identical(
    x$dist,
    "poisson"
  )
})


test_that("frequency specifications can carry labels", {
  x <- frequency_spec(
    "poisson",
    lambda = 8,
    label = "Annual motor claims"
  )

  expect_identical(
    x$label,
    "Annual motor claims"
  )
})


test_that("frequency specifications do not store computational state", {
  x <- frequency_spec(
    "poisson",
    lambda = 5
  )

  forbidden <- c(
    "seed",
    "n_sim",
    "rng",
    "pmf",
    "cdf",
    "quantile",
    "raw_moment",
    "factorial_moment",
    "pgf",
    "panjer"
  )

  expect_false(
    any(
      forbidden %in% names(x)
    )
  )
})


# dist validation ----------------------------------------------------------

test_that("dist must be a single non-missing string", {
  expect_error(
    frequency_spec(
      NULL,
      lambda = 1
    ),
    "`dist`"
  )

  expect_error(
    frequency_spec(
      character(),
      lambda = 1
    ),
    "`dist`"
  )

  expect_error(
    frequency_spec(
      c("poisson", "binomial"),
      lambda = 1
    ),
    "`dist`"
  )

  expect_error(
    frequency_spec(
      NA_character_,
      lambda = 1
    ),
    "`dist`"
  )

  expect_error(
    frequency_spec(
      "",
      lambda = 1
    ),
    "`dist`"
  )

  expect_error(
    frequency_spec(
      "   ",
      lambda = 1
    ),
    "`dist`"
  )

  expect_error(
    frequency_spec(
      1,
      lambda = 1
    ),
    "`dist`"
  )
})


test_that("unsupported distributions are rejected", {
  expect_error(
    frequency_spec(
      "normal",
      mean = 0,
      sd = 1
    ),
    "must be one of"
  )
})


# Parameter naming ---------------------------------------------------------

test_that("all distribution parameters must be named", {
  expect_error(
    frequency_spec(
      "poisson",
      4
    ),
    "must be named"
  )
})


test_that("duplicated parameter names are rejected", {
  expect_error(
    frequency_spec(
      "poisson",
      lambda = 2,
      lambda = 3
    ),
    "must be unique"
  )
})


test_that("missing parameters are rejected", {
  expect_error(
    frequency_spec("poisson"),
    "Missing required parameter"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = 10
    ),
    "Missing required parameter"
  )

  expect_error(
    frequency_spec(
      "negative_binomial",
      prob = 0.4
    ),
    "Missing required parameter"
  )

  expect_error(
    frequency_spec("geometric"),
    "Missing required parameter"
  )
})


test_that("unused parameters are rejected", {
  expect_error(
    frequency_spec(
      "poisson",
      lambda = 4,
      prob = 0.5
    ),
    "Unused parameter"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = 10,
      prob = 0.2,
      lambda = 4
    ),
    "Unused parameter"
  )
})


# Label validation ---------------------------------------------------------

test_that("label validation is strict", {
  expect_error(
    frequency_spec(
      "poisson",
      lambda = 1,
      label = ""
    ),
    "`label`"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = 1,
      label = "   "
    ),
    "`label`"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = 1,
      label = NA_character_
    ),
    "`label`"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = 1,
      label = c("a", "b")
    ),
    "`label`"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = 1,
      label = 123
    ),
    "`label`"
  )
})


# Poisson validation -------------------------------------------------------

test_that("Poisson accepts its full valid boundary", {
  expect_no_error(
    frequency_spec(
      "poisson",
      lambda = 0
    )
  )

  expect_no_error(
    frequency_spec(
      "poisson",
      lambda = 1e-12
    )
  )

  expect_no_error(
    frequency_spec(
      "poisson",
      lambda = 1e6
    )
  )
})


test_that("Poisson rejects invalid lambda", {
  expect_error(
    frequency_spec(
      "poisson",
      lambda = -1
    ),
    "lambda"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = Inf
    ),
    "finite numeric scalar"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = NA_real_
    ),
    "finite numeric scalar"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = c(1, 2)
    ),
    "finite numeric scalar"
  )

  expect_error(
    frequency_spec(
      "poisson",
      lambda = "4"
    ),
    "finite numeric scalar"
  )
})


# Binomial validation ------------------------------------------------------

test_that("binomial accepts valid boundary cases", {
  expect_no_error(
    frequency_spec(
      "binomial",
      size = 0,
      prob = 0.5
    )
  )

  expect_no_error(
    frequency_spec(
      "binomial",
      size = 10,
      prob = 0
    )
  )

  expect_no_error(
    frequency_spec(
      "binomial",
      size = 10,
      prob = 1
    )
  )
})


test_that("binomial rejects invalid size", {
  expect_error(
    frequency_spec(
      "binomial",
      size = -1,
      prob = 0.5
    ),
    "size"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = 2.5,
      prob = 0.5
    ),
    "integer-valued"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = Inf,
      prob = 0.5
    ),
    "finite numeric scalar"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = NA_real_,
      prob = 0.5
    ),
    "finite numeric scalar"
  )
})


test_that("binomial rejects invalid probability", {
  expect_error(
    frequency_spec(
      "binomial",
      size = 10,
      prob = -0.01
    ),
    "prob"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = 10,
      prob = 1.01
    ),
    "prob"
  )

  expect_error(
    frequency_spec(
      "binomial",
      size = 10,
      prob = NA_real_
    ),
    "finite numeric scalar"
  )
})


# Negative binomial validation --------------------------------------------

test_that("negative binomial accepts non-integer positive size", {
  expect_no_error(
    frequency_spec(
      "negative_binomial",
      size = 0.1,
      prob = 0.5
    )
  )

  expect_no_error(
    frequency_spec(
      "negative_binomial",
      size = 2.5,
      prob = 0.5
    )
  )

  expect_no_error(
    frequency_spec(
      "negative_binomial",
      size = 5,
      prob = 1
    )
  )
})


test_that("negative binomial rejects invalid size", {
  expect_error(
    frequency_spec(
      "negative_binomial",
      size = 0,
      prob = 0.5
    ),
    "size"
  )

  expect_error(
    frequency_spec(
      "negative_binomial",
      size = -1,
      prob = 0.5
    ),
    "size"
  )

  expect_error(
    frequency_spec(
      "negative_binomial",
      size = Inf,
      prob = 0.5
    ),
    "finite numeric scalar"
  )
})


test_that("negative binomial requires probability in (0, 1]", {
  expect_error(
    frequency_spec(
      "negative_binomial",
      size = 2,
      prob = 0
    ),
    "prob"
  )

  expect_error(
    frequency_spec(
      "negative_binomial",
      size = 2,
      prob = -0.1
    ),
    "prob"
  )

  expect_error(
    frequency_spec(
      "negative_binomial",
      size = 2,
      prob = 1.1
    ),
    "prob"
  )
})


# Geometric validation -----------------------------------------------------

test_that("geometric accepts probability in (0, 1]", {
  expect_no_error(
    frequency_spec(
      "geometric",
      prob = 1e-10
    )
  )

  expect_no_error(
    frequency_spec(
      "geometric",
      prob = 1
    )
  )
})


test_that("geometric rejects invalid probability", {
  expect_error(
    frequency_spec(
      "geometric",
      prob = 0
    ),
    "prob"
  )

  expect_error(
    frequency_spec(
      "geometric",
      prob = -0.1
    ),
    "prob"
  )

  expect_error(
    frequency_spec(
      "geometric",
      prob = 1.1
    ),
    "prob"
  )
})


# Registry ----------------------------------------------------------------

test_that("built-in frequency distributions expose expected capabilities", {
  registry <- frequency_registry()

  expected <- c(
    "validate",
    "support",
    "pmf",
    "cdf",
    "quantile",
    "rng",
    "raw_moment",
    "factorial_moment",
    "pgf",
    "panjer"
  )

  for (
    dist in c(
      "poisson",
      "binomial",
      "negative_binomial",
      "geometric"
    )
  ) {
    expect_true(
      all(
        expected %in%
          names(registry[[dist]])
      )
    )

    for (operation in expected) {
      expect_true(
        is.function(
    registry[[dist]][[operation]]
        )
      )
    }
  }
})


# PMF ---------------------------------------------------------------------

test_that("Poisson pmf agrees with stats", {
  entry <- frequency_registry()$poisson
  params <- list(lambda = 3.7)
  x <- 0:12

  expect_equal(
    entry$pmf(x, params),
    stats::dpois(
      x,
      lambda = 3.7
    )
  )

  expect_equal(
    entry$pmf(
      x,
      params,
      log = TRUE
    ),
    stats::dpois(
      x,
      lambda = 3.7,
      log = TRUE
    )
  )
})


test_that("binomial pmf agrees with stats", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 12,
    prob = 0.35
  )

  x <- 0:12

  expect_equal(
    entry$pmf(x, params),
    stats::dbinom(
      x,
      size = 12,
      prob = 0.35
    )
  )
})


test_that("negative binomial pmf agrees with stats", {
  entry <- frequency_registry()$negative_binomial

  params <- list(
    size = 2.5,
    prob = 0.4
  )

  x <- 0:20

  expect_equal(
    entry$pmf(x, params),
    stats::dnbinom(
      x,
      size = 2.5,
      prob = 0.4
    )
  )
})


test_that("geometric pmf agrees with stats", {
  entry <- frequency_registry()$geometric

  params <- list(
    prob = 0.3
  )

  x <- 0:20

  expect_equal(
    entry$pmf(x, params),
    stats::dgeom(
      x,
      prob = 0.3
    )
  )
})


# PMF normalization --------------------------------------------------------

test_that("binomial pmf sums exactly to one numerically", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 40,
    prob = 0.18
  )

  expect_equal(
    sum(
      entry$pmf(
        0:40,
        params
      )
    ),
    1,
    tolerance = 1e-14
  )
})


test_that("infinite-support pmfs have total probability one", {
  registry <- frequency_registry()

  poisson <- registry$poisson
  poisson_params <- list(
    lambda = 4
  )

  expect_equal(
    sum(
      poisson$pmf(
        0:100,
        poisson_params
      )
    ),
    1,
    tolerance = 1e-12
  )

  nb <- registry$negative_binomial
  nb_params <- list(
    size = 3,
    prob = 0.6
  )

  expect_equal(
    sum(
      nb$pmf(
        0:500,
        nb_params
      )
    ),
    1,
    tolerance = 1e-12
  )

  geom <- registry$geometric
  geom_params <- list(
    prob = 0.4
  )

  expect_equal(
    sum(
      geom$pmf(
        0:500,
        geom_params
      )
    ),
    1,
    tolerance = 1e-12
  )
})


# CDF ---------------------------------------------------------------------

test_that("built-in CDFs agree with stats", {
  registry <- frequency_registry()

  q <- 0:10

  expect_equal(
    registry$poisson$cdf(
      q,
      list(lambda = 4)
    ),
    stats::ppois(
      q,
      lambda = 4
    )
  )

  expect_equal(
    registry$binomial$cdf(
      q,
      list(
        size = 12,
        prob = 0.3
      )
    ),
    stats::pbinom(
      q,
      size = 12,
      prob = 0.3
    )
  )

  expect_equal(
    registry$negative_binomial$cdf(
      q,
      list(
        size = 2.5,
        prob = 0.45
      )
    ),
    stats::pnbinom(
      q,
      size = 2.5,
      prob = 0.45
    )
  )

  expect_equal(
    registry$geometric$cdf(
      q,
      list(prob = 0.4)
    ),
    stats::pgeom(
      q,
      prob = 0.4
    )
  )
})


test_that("CDF tail and log options are preserved", {
  entry <- frequency_registry()$poisson

  params <- list(
    lambda = 4
  )

  q <- 0:10

  expect_equal(
    entry$cdf(
      q,
      params,
      lower.tail = FALSE,
      log.p = TRUE
    ),
    stats::ppois(
      q,
      lambda = 4,
      lower.tail = FALSE,
      log.p = TRUE
    )
  )
})


# Quantiles ---------------------------------------------------------------

test_that("built-in quantiles agree with stats", {
  registry <- frequency_registry()

  p <- c(
    0,
    0.1,
    0.25,
    0.5,
    0.9,
    0.99,
    1
  )

  expect_equal(
    registry$poisson$quantile(
      p,
      list(lambda = 4)
    ),
    stats::qpois(
      p,
      lambda = 4
    )
  )

  expect_equal(
    registry$binomial$quantile(
      p,
      list(
        size = 20,
        prob = 0.25
      )
    ),
    stats::qbinom(
      p,
      size = 20,
      prob = 0.25
    )
  )

  expect_equal(
    registry$negative_binomial$quantile(
      p,
      list(
        size = 3,
        prob = 0.5
      )
    ),
    stats::qnbinom(
      p,
      size = 3,
      prob = 0.5
    )
  )

  expect_equal(
    registry$geometric$quantile(
      p,
      list(prob = 0.4)
    ),
    stats::qgeom(
      p,
      prob = 0.4
    )
  )
})


# RNG ---------------------------------------------------------------------

test_that("Poisson RNG delegates consistently to stats", {
  entry <- frequency_registry()$poisson
  params <- list(lambda = 4)

  set.seed(1234)
  obtained <- entry$rng(
    20,
    params
  )

  set.seed(1234)
  expected <- stats::rpois(
    20,
    lambda = 4
  )

  expect_identical(
    obtained,
    expected
  )
})


test_that("binomial RNG delegates consistently to stats", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 20,
    prob = 0.3
  )

  set.seed(1234)
  obtained <- entry$rng(
    20,
    params
  )

  set.seed(1234)
  expected <- stats::rbinom(
    20,
    size = 20,
    prob = 0.3
  )

  expect_identical(
    obtained,
    expected
  )
})


test_that("negative binomial RNG delegates consistently to stats", {
  entry <- frequency_registry()$negative_binomial

  params <- list(
    size = 2.5,
    prob = 0.4
  )

  set.seed(1234)
  obtained <- entry$rng(
    20,
    params
  )

  set.seed(1234)
  expected <- stats::rnbinom(
    20,
    size = 2.5,
    prob = 0.4
  )

  expect_identical(
    obtained,
    expected
  )
})


test_that("geometric RNG delegates consistently to stats", {
  entry <- frequency_registry()$geometric

  params <- list(
    prob = 0.4
  )

  set.seed(1234)
  obtained <- entry$rng(
    20,
    params
  )

  set.seed(1234)
  expected <- stats::rgeom(
    20,
    prob = 0.4
  )

  expect_identical(
    obtained,
    expected
  )
})


# Factorial moments --------------------------------------------------------

test_that("Poisson factorial moments are correct", {
  entry <- frequency_registry()$poisson

  params <- list(
    lambda = 4
  )

  for (order in 0:6) {
    expect_equal(
      entry$factorial_moment(
        order,
        params
      ),
      4^order,
      tolerance = 1e-12
    )
  }
})


test_that("binomial factorial moments are correct", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 10,
    prob = 0.3
  )

  expected <- c(
    1,
    10 * 0.3,
    10 * 9 * 0.3^2,
    10 * 9 * 8 * 0.3^3,
    10 * 9 * 8 * 7 * 0.3^4
  )

  for (order in 0:4) {
    expect_equal(
      entry$factorial_moment(
        order,
        params
      ),
      expected[order + 1L],
      tolerance = 1e-12
    )
  }

  expect_equal(
    entry$factorial_moment(
      11,
      params
    ),
    0
  )
})


test_that("negative binomial factorial moments are correct", {
  entry <- frequency_registry()$negative_binomial

  r <- 2.5
  p <- 0.4
  q <- 1 - p

  params <- list(
    size = r,
    prob = p
  )

  expected <- c(
    1,
    r * q / p,
    r * (r + 1) * (q / p)^2,
    r * (r + 1) * (r + 2) * (q / p)^3,
    r * (r + 1) * (r + 2) * (r + 3) * (q / p)^4
  )

  for (order in 0:4) {
    expect_equal(
      entry$factorial_moment(
        order,
        params
      ),
      expected[order + 1L],
      tolerance = 1e-12
    )
  }
})


test_that("geometric factorial moments are correct", {
  entry <- frequency_registry()$geometric

  p <- 0.4
  q <- 1 - p

  params <- list(
    prob = p
  )

  for (order in 0:6) {
    expected <- factorial(order) *
      (q / p)^order

    expect_equal(
      entry$factorial_moment(
        order,
        params
      ),
      expected,
      tolerance = 1e-12
    )
  }
})


# Raw moments --------------------------------------------------------------

test_that("Poisson raw moments match known formulas", {
  entry <- frequency_registry()$poisson

  lambda <- 4

  params <- list(
    lambda = lambda
  )

  expected <- c(
    1,
    lambda,
    lambda^2 + lambda,
    lambda^3 + 3 * lambda^2 + lambda,
    lambda^4 +
      6 * lambda^3 +
      7 * lambda^2 +
      lambda
  )

  for (order in 0:4) {
    expect_equal(
      entry$raw_moment(
        order,
        params
      ),
      expected[order + 1L],
      tolerance = 1e-12
    )
  }
})


test_that("binomial raw moments agree with direct summation", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 15,
    prob = 0.27
  )

  x <- 0:15
  probabilities <- stats::dbinom(
    x,
    size = 15,
    prob = 0.27
  )

  for (order in 0:4) {
    expected <- sum(
      x^order *
        probabilities
    )

    expect_equal(
      entry$raw_moment(
        order,
        params
      ),
      expected,
      tolerance = 1e-11
    )
  }
})


test_that("negative binomial raw moments agree with direct summation", {
  entry <- frequency_registry()$negative_binomial

  params <- list(
    size = 3.5,
    prob = 0.7
  )

  x <- 0:500

  probabilities <- stats::dnbinom(
    x,
    size = params$size,
    prob = params$prob
  )

  for (order in 0:4) {
    expected <- sum(
      x^order *
        probabilities
    )

    expect_equal(
      entry$raw_moment(
        order,
        params
      ),
      expected,
      tolerance = 1e-9
    )
  }
})


test_that("geometric raw moments agree with direct summation", {
  entry <- frequency_registry()$geometric

  params <- list(
    prob = 0.4
  )

  x <- 0:500

  probabilities <- stats::dgeom(
    x,
    prob = params$prob
  )

  for (order in 0:4) {
    expected <- sum(
      x^order *
        probabilities
    )

    expect_equal(
      entry$raw_moment(
        order,
        params
      ),
      expected,
      tolerance = 1e-9
    )
  }
})


test_that("first and second moments recover known means and variances", {
  registry <- frequency_registry()

  # Poisson
  lambda <- 5.3
  params <- list(lambda = lambda)

  m1 <- registry$poisson$raw_moment(
    1,
    params
  )

  m2 <- registry$poisson$raw_moment(
    2,
    params
  )

  expect_equal(
    m1,
    lambda
  )

  expect_equal(
    m2 - m1^2,
    lambda,
    tolerance = 1e-12
  )

  # Binomial
  m <- 20
  p <- 0.3

  params <- list(
    size = m,
    prob = p
  )

  m1 <- registry$binomial$raw_moment(
    1,
    params
  )

  m2 <- registry$binomial$raw_moment(
    2,
    params
  )

  expect_equal(
    m1,
    m * p,
    tolerance = 1e-12
  )

  expect_equal(
    m2 - m1^2,
    m * p * (1 - p),
    tolerance = 1e-12
  )

  # Negative binomial
  r <- 2.5
  p <- 0.4

  params <- list(
    size = r,
    prob = p
  )

  m1 <- registry$negative_binomial$raw_moment(
    1,
    params
  )

  m2 <- registry$negative_binomial$raw_moment(
    2,
    params
  )

  expect_equal(
    m1,
    r * (1 - p) / p,
    tolerance = 1e-12
  )

  expect_equal(
    m2 - m1^2,
    r * (1 - p) / p^2,
    tolerance = 1e-12
  )

  # Geometric
  p <- 0.35

  params <- list(
    prob = p
  )

  m1 <- registry$geometric$raw_moment(
    1,
    params
  )

  m2 <- registry$geometric$raw_moment(
    2,
    params
  )

  expect_equal(
    m1,
    (1 - p) / p,
    tolerance = 1e-12
  )

  expect_equal(
    m2 - m1^2,
    (1 - p) / p^2,
    tolerance = 1e-12
  )
})


# Moment order validation --------------------------------------------------

test_that("moment orders must be nonnegative integers", {
  entry <- frequency_registry()$poisson

  params <- list(lambda = 2)

  expect_error(
    entry$raw_moment(
      -1,
      params
    ),
    "nonnegative integer"
  )

  expect_error(
    entry$raw_moment(
      1.5,
      params
    ),
    "nonnegative integer"
  )

  expect_error(
    entry$raw_moment(
      NA_real_,
      params
    ),
    "nonnegative integer"
  )

  expect_error(
    entry$raw_moment(
      Inf,
      params
    ),
    "nonnegative integer"
  )

  expect_error(
    entry$raw_moment(
      c(1, 2),
      params
    ),
    "nonnegative integer"
  )
})


# Stirling numbers ---------------------------------------------------------

test_that("internal Stirling numbers of the second kind are correct", {
  S2 <- tidyactuarial:::.frequency_stirling2

  expect_equal(
    S2(0, 0),
    1
  )

  expect_equal(
    S2(1, 1),
    1
  )

  expect_equal(
    S2(2, 1),
    1
  )

  expect_equal(
    S2(2, 2),
    1
  )

  expect_equal(
    S2(3, 1),
    1
  )

  expect_equal(
    S2(3, 2),
    3
  )

  expect_equal(
    S2(3, 3),
    1
  )

  expect_equal(
    S2(4, 1),
    1
  )

  expect_equal(
    S2(4, 2),
    7
  )

  expect_equal(
    S2(4, 3),
    6
  )

  expect_equal(
    S2(4, 4),
    1
  )

  expect_equal(
    S2(3, 4),
    0
  )
})


# PGFs --------------------------------------------------------------------

test_that("all built-in PGFs satisfy G(1) = 1", {
  registry <- frequency_registry()

  cases <- list(
    poisson = list(lambda = 3.2),
    binomial = list(
      size = 12,
      prob = 0.3
    ),
    negative_binomial = list(
      size = 2.5,
      prob = 0.4
    ),
    geometric = list(
      prob = 0.4
    )
  )

  for (dist in names(cases)) {
    expect_equal(
      registry[[dist]]$pgf(
        1,
        cases[[dist]]
      ),
      1,
      tolerance = 1e-12
    )
  }
})


test_that("PGFs agree with sums of z^N p_N", {
  registry <- frequency_registry()

  z <- 0.7

  # Poisson
  params <- list(lambda = 4)
  x <- 0:100

  expected <- sum(
    z^x *
      registry$poisson$pmf(
        x,
        params
      )
  )

  expect_equal(
    registry$poisson$pgf(
      z,
      params
    ),
    expected,
    tolerance = 1e-12
  )

  # Binomial
  params <- list(
    size = 20,
    prob = 0.25
  )

  x <- 0:20

  expected <- sum(
    z^x *
      registry$binomial$pmf(
        x,
        params
      )
  )

  expect_equal(
    registry$binomial$pgf(
      z,
      params
    ),
    expected,
    tolerance = 1e-12
  )

  # Negative binomial
  params <- list(
    size = 3,
    prob = 0.6
  )

  x <- 0:500

  expected <- sum(
    z^x *
      registry$negative_binomial$pmf(
        x,
        params
      )
  )

  expect_equal(
    registry$negative_binomial$pgf(
      z,
      params
    ),
    expected,
    tolerance = 1e-12
  )

  # Geometric
  params <- list(
    prob = 0.4
  )

  x <- 0:500

  expected <- sum(
    z^x *
      registry$geometric$pmf(
        x,
        params
      )
  )

  expect_equal(
    registry$geometric$pgf(
      z,
      params
    ),
    expected,
    tolerance = 1e-12
  )
})


# Panjer ------------------------------------------------------------------

test_that("Poisson satisfies the Panjer recursion", {
  expect_panjer_recursion(
    entry = frequency_registry()$poisson,
    params = list(lambda = 4),
    k = 1:20
  )
})


test_that("binomial satisfies the Panjer recursion", {
  expect_panjer_recursion(
    entry = frequency_registry()$binomial,
    params = list(
      size = 12,
      prob = 0.3
    ),
    k = 1:12
  )
})


test_that("negative binomial satisfies the Panjer recursion", {
  expect_panjer_recursion(
    entry = frequency_registry()$negative_binomial,
    params = list(
      size = 2.5,
      prob = 0.4
    ),
    k = 1:20
  )
})


test_that("geometric satisfies the Panjer recursion", {
  expect_panjer_recursion(
    entry = frequency_registry()$geometric,
    params = list(
      prob = 0.4
    ),
    k = 1:20
  )
})


test_that("binomial prob = 1 remains valid but has no regular Panjer representation", {
  x <- frequency_spec(
    "binomial",
    size = 10,
    prob = 1
  )

  expect_valid_frequency_spec(
    x,
    "binomial"
  )

  panjer <- frequency_registry()$binomial$panjer(
    x$params
  )

  expect_null(panjer)
})


# Degenerate boundary distributions ---------------------------------------

test_that("Poisson lambda = 0 is degenerate at zero", {
  entry <- frequency_registry()$poisson

  params <- list(
    lambda = 0
  )

  expect_equal(
    entry$pmf(
      0,
      params
    ),
    1
  )

  expect_equal(
    entry$pmf(
      1:10,
      params
    ),
    rep(
      0,
      10
    )
  )

  expect_equal(
    entry$raw_moment(
      1,
      params
    ),
    0
  )

  expect_equal(
    entry$factorial_moment(
      4,
      params
    ),
    0
  )
})


test_that("binomial size = 0 is degenerate at zero", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 0,
    prob = 0.7
  )

  expect_equal(
    entry$pmf(
      0,
      params
    ),
    1
  )

  expect_equal(
    entry$raw_moment(
      1,
      params
    ),
    0
  )
})


test_that("binomial prob = 0 is degenerate at zero", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 10,
    prob = 0
  )

  expect_equal(
    entry$pmf(
      0,
      params
    ),
    1
  )

  expect_equal(
    entry$pmf(
      1:10,
      params
    ),
    rep(
      0,
      10
    )
  )

  expect_equal(
    entry$raw_moment(
      1,
      params
    ),
    0
  )
})


test_that("binomial prob = 1 is degenerate at size", {
  entry <- frequency_registry()$binomial

  params <- list(
    size = 10,
    prob = 1
  )

  expect_equal(
    entry$pmf(
      10,
      params
    ),
    1
  )

  expect_equal(
    entry$pmf(
      0:9,
      params
    ),
    rep(
      0,
      10
    )
  )

  expect_equal(
    entry$raw_moment(
      1,
      params
    ),
    10
  )
})


test_that("negative binomial prob = 1 is degenerate at zero", {
  entry <- frequency_registry()$negative_binomial

  params <- list(
    size = 2.5,
    prob = 1
  )

  expect_equal(
    entry$pmf(
      0,
      params
    ),
    1
  )

  expect_equal(
    entry$raw_moment(
      1,
      params
    ),
    0
  )
})


test_that("geometric prob = 1 is degenerate at zero", {
  entry <- frequency_registry()$geometric

  params <- list(
    prob = 1
  )

  expect_equal(
    entry$pmf(
      0,
      params
    ),
    1
  )

  expect_equal(
    entry$raw_moment(
      1,
      params
    ),
    0
  )
})


# Custom distributions -----------------------------------------------------

test_that("custom frequency specification requires a named list", {
  expect_error(
    frequency_spec(
      "custom",
      custom = NULL
    ),
    "named list"
  )

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        function(x, params) 1
      )
    ),
    "named list"
  )
})


test_that("custom frequency specification requires a pmf", {
  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        cdf = function(q, params) {
          rep(
            1,
            length(q)
          )
        }
      )
    ),
    "custom\\$pmf"
  )
})


test_that("custom pmf must be a function", {
  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = 1
      )
    ),
    "must be supplied as a function"
  )
})


test_that("optional custom operations must be functions", {
  pmf <- function(x, params, log = FALSE) {
    rep(
      0.5,
      length(x)
    )
  }

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        cdf = 1
      )
    ),
    "custom\\$cdf"
  )

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        rng = "not a function"
      )
    ),
    "custom\\$rng"
  )

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        pgf = 1
      )
    ),
    "custom\\$pgf"
  )
})


test_that("unknown custom components are rejected", {
  pmf <- function(x, params, log = FALSE) {
    rep(
      0.5,
      length(x)
    )
  }

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        something_else = function(...) 1
      )
    ),
    "Unsupported component"
  )
})


test_that("custom frequency uses default count support", {
  pmf <- function(x, params, log = FALSE) {
    p <- ifelse(
      x == 0,
      1,
      0
    )

    if (log) {
      log(p)
    } else {
      p
    }
  }

  x <- frequency_spec(
    "custom",
    custom = list(
      pmf = pmf
    ),
    label = "Degenerate custom frequency"
  )

  expect_valid_frequency_spec(
    x,
    "custom"
  )

  expect_count_support(
    x,
    lower = 0,
    upper = Inf,
    step = 1
  )

  expect_true(
    is.function(
    x$custom$pmf
    )
  )
})


test_that("custom frequency can use finite support", {
  pmf <- function(x, params, log = FALSE) {
    probs <- params$probabilities
    support <- seq_along(probs) - 1L

    out <- numeric(
      length(x)
    )

    index <- match(
      x,
      support
    )

    valid <- !is.na(index)

    out[valid] <- probs[
      index[valid]
    ]

    if (log) {
      log(out)
    } else {
      out
    }
  }

  x <- frequency_spec(
    "custom",
    probabilities = c(
      0.5,
      0.3,
      0.2
    ),
    custom = list(
      pmf = pmf,
      support = list(
        lower = 0,
        upper = 2,
        step = 1
      )
    )
  )

  expect_valid_frequency_spec(
    x,
    "custom"
  )

  expect_count_support(
    x,
    lower = 0,
    upper = 2,
    step = 1
  )

  expect_equal(
    x$params$probabilities,
    c(
      0.5,
      0.3,
      0.2
    )
  )

  expect_null(
    x$custom$support
  )
})


test_that("custom support may depend on parameters", {
  pmf <- function(x, params, log = FALSE) {
    stats::dbinom(
      x,
      size = params$maximum,
      prob = params$prob,
      log = log
    )
  }

  support <- function(params) {
    list(
      type = "count",
      lower = 0,
      upper = params$maximum,
      step = 1
    )
  }

  x <- frequency_spec(
    "custom",
    maximum = 7,
    prob = 0.4,
    custom = list(
      pmf = pmf,
      support = support
    )
  )

  expect_count_support(
    x,
    lower = 0,
    upper = 7,
    step = 1
  )
})


test_that("custom support must describe a valid count grid", {
  pmf <- function(x, params, log = FALSE) {
    rep(
      0,
      length(x)
    )
  }

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        support = list(
          lower = -1,
          upper = 5,
          step = 1
        )
      )
    ),
    "support\\$lower"
  )

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        support = list(
          lower = 0,
          upper = 5,
          step = 0
        )
      )
    ),
    "support\\$step"
  )

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        support = list(
          lower = 5,
          upper = 3,
          step = 1
        )
      )
    ),
    "support\\$upper"
  )

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        support = list(
          type = "continuous",
          lower = 0,
          upper = 5,
          step = 1
        )
      )
    ),
    "type ="
  )
})


test_that("custom support upper endpoint must belong to the grid", {
  pmf <- function(x, params, log = FALSE) {
    rep(
      0,
      length(x)
    )
  }

  expect_error(
    frequency_spec(
      "custom",
      custom = list(
        pmf = pmf,
        support = list(
          lower = 0,
          upper = 5,
          step = 2
        )
      )
    ),
    "support grid"
  )
})


test_that("custom cannot be supplied for built-in distributions", {
  expect_error(
    frequency_spec(
      "poisson",
      lambda = 4,
      custom = list(
        pmf = function(...) 1
      )
    ),
    "can only be supplied"
  )
})


# Print method -------------------------------------------------------------

test_that("print method returns the object invisibly", {
  x <- frequency_spec(
    "poisson",
    lambda = 4
  )

  result <- invisible(
    capture.output(
      returned <- print(x)
    )
  )

  expect_identical(
    returned,
    x
  )
})


test_that("print method displays core specification information", {
  x <- frequency_spec(
    "poisson",
    lambda = 4,
    label = "Annual claims"
  )

  expect_output(
    print(x),
    "<tidyact_frequency_spec>",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Distribution: Poisson",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "lambda: 4",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Support:",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Label: Annual claims",
    fixed = TRUE
  )
})


test_that("finite support is represented in print output", {
  x <- frequency_spec(
    "binomial",
    size = 10,
    prob = 0.3
  )

  expect_output(
    print(x),
    "Support: {0, 1, ..., 10}",
    fixed = TRUE
  )
})
