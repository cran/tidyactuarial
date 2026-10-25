# Tests for discretize_severity ---------------------------------------------

.disc_full_grid <- function(x) {
  step <- x$support$step
  to <- x$support$upper

  grid <- seq(
    0,
    to,
    by = step
  )

  probs <- numeric(
    length(grid)
  )

  idx <- match(
    x$params$values,
    grid
  )

  probs[idx] <- x$params$probs

  list(
    grid = grid,
    probs = probs
  )
}


test_that("discretize_severity returns a lattice severity specification", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1,
    label = "Original severity"
  )

  out <- discretize_severity(
    x,
    step = 0.5,
    to = 6,
    method = "upper"
  )

  expect_s3_class(
    out,
    "tidyact_severity_spec"
  )

  expect_s3_class(
    out,
    "tidyact_distribution_spec"
  )

  expect_identical(
    out$schema_version,
    1L
  )

  expect_identical(
    out$dist,
    "discrete"
  )

  expect_identical(
    out$support$type,
    "lattice"
  )

  expect_equal(
    out$support$lower,
    0
  )

  expect_equal(
    out$support$upper,
    6
  )

  expect_equal(
    out$support$step,
    0.5
  )

  expect_identical(
    out$label,
    "Original severity"
  )

  expect_null(out$data)
  expect_null(out$custom)
})


test_that("all discretization methods produce proper probabilities", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  methods <- c(
    "upper",
    "lower",
    "rounding",
    "unbiased"
  )

  for (method in methods) {
    out <- discretize_severity(
      x,
      step = 0.5,
      to = 8,
      method = method
    )

    expect_true(
      all(
        is.finite(
          out$params$probs
        )
      )
    )

    expect_true(
      all(
        out$params$probs > 0
      )
    )

    expect_equal(
      sum(
        out$params$probs
      ),
      1,
      tolerance = 1e-12
    )

    expect_true(
      all(
        out$params$values >= 0
      )
    )

    expect_true(
      all(
        out$params$values <= 8
      )
    )
  }
})


test_that("upper discretization follows forward CDF differences", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = 1,
    to = 3,
    method = "upper"
  )

  full <- .disc_full_grid(out)

  F <- function(z) {
    stats::pexp(
      z,
      rate = 1
    )
  }

  expected <- c(
    F(1) - F(0),
    F(2) - F(1),
    F(3) - F(2),
    1 - F(3)
  )

  expect_equal(
    full$grid,
    0:3
  )

  expect_equal(
    full$probs,
    expected,
    tolerance = 1e-13
  )
})


test_that("lower discretization follows backward CDF differences", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = 1,
    to = 3,
    method = "lower"
  )

  full <- .disc_full_grid(out)

  F <- function(z) {
    stats::pexp(
      z,
      rate = 1
    )
  }

  expected <- c(
    F(0),
    F(1) - F(0),
    F(2) - F(1),
    F(3) - F(2) +
      1 - F(3)
  )

  expect_equal(
    full$probs,
    expected,
    tolerance = 1e-13
  )

  expect_equal(
    full$probs[4],
    1 - F(2),
    tolerance = 1e-13
  )
})


test_that("rounding discretization follows midpoint cells", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = 1,
    to = 3,
    method = "rounding"
  )

  full <- .disc_full_grid(out)

  F <- function(z) {
    stats::pexp(
      z,
      rate = 1
    )
  }

  expected <- c(
    F(0.5),
    F(1.5) - F(0.5),
    F(2.5) - F(1.5),
    1 - F(2.5)
  )

  expect_equal(
    full$probs,
    expected,
    tolerance = 1e-13
  )
})


test_that("unbiased discretization matches the limited first moment", {
  rate <- 1

  x <- severity_spec(
    "exponential",
    rate = rate
  )

  out <- discretize_severity(
    x,
    step = 0.5,
    to = 4,
    method = "unbiased"
  )

  full <- .disc_full_grid(out)

  discretized_mean <- sum(
    full$grid *
      full$probs
  )

  limited_mean <- (
    1 -
      exp(
        -rate * 4
      )
  ) /
    rate

  expect_equal(
    discretized_mean,
    limited_mean,
    tolerance = 1e-11
  )

  expect_identical(
    out$metadata$discretization$limited_mean_source,
    "analytic"
  )
})


test_that("unbiased discretization agrees with its analytical cell formulas", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  h <- 1
  b <- 3

  out <- discretize_severity(
    x,
    step = h,
    to = b,
    method = "unbiased"
  )

  full <- .disc_full_grid(out)

  L <- function(d) {
    1 -
      exp(
        -d
      )
  }

  expected <- c(
    1 - L(1),
    2 * L(1) -
      L(0) -
      L(2),
    2 * L(2) -
      L(1) -
      L(3),
    L(3) -
      L(2)
  )

  expect_equal(
    full$probs,
    expected,
    tolerance = 1e-12
  )

  expect_equal(
    sum(expected),
    1,
    tolerance = 1e-12
  )
})


test_that("automatic cutoff is a lattice point with controlled tail probability", {
  tail_prob <- 1e-5
  step <- 0.25

  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = step,
    tail_prob = tail_prob,
    method = "upper"
  )

  info <- out$metadata$discretization

  expect_identical(
    info$cutoff_source,
    "tail_prob"
  )

  expect_equal(
    info$requested_tail_prob,
    tail_prob
  )

  expect_lte(
    info$actual_tail_probability,
    tail_prob +
      1e-12
  )

  ratio <- info$to /
    step

  expect_equal(
    ratio,
    round(ratio),
    tolerance = 1e-12
  )
})


test_that("supplied cutoff is recorded explicitly", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = 0.5,
    to = 5,
    method = "rounding"
  )

  info <- out$metadata$discretization

  expect_identical(
    info$cutoff_source,
    "supplied"
  )

  expect_true(
    is.na(
      info$requested_tail_prob
    )
  )

  expect_equal(
    info$to,
    5
  )

  expect_identical(
    info$tail_handling,
    "assigned_to_upper_grid_point"
  )
})


test_that("tail mass is retained rather than silently discarded", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = 1,
    to = 2,
    method = "lower"
  )

  full <- .disc_full_grid(out)

  tail <- exp(-2)

  expect_equal(
    out$metadata$discretization$actual_tail_probability,
    tail,
    tolerance = 1e-13
  )

  expect_equal(
    sum(full$probs),
    1,
    tolerance = 1e-13
  )

  expect_equal(
    full$probs[3],
    exp(-1),
    tolerance = 1e-13
  )
})


test_that("discretization metadata preserves the source distribution", {
  x <- severity_spec(
    "lognormal",
    meanlog = 1,
    sdlog = 0.5
  )

  out <- discretize_severity(
    x,
    step = 0.5,
    to = 10,
    method = "unbiased",
    label = "Discretized loss"
  )

  expect_identical(
    out$metadata$discretization$source_dist,
    "lognormal"
  )

  expect_identical(
    out$metadata$discretization$method,
    "unbiased"
  )

  expect_equal(
    out$metadata$discretization$step,
    0.5
  )

  expect_identical(
    out$label,
    "Discretized loss"
  )

  expect_null(
    x$metadata$discretization
  )
})


test_that("Burr unbiased discretization records numerical limited-mean evaluation", {
  x <- severity_spec(
    "burr",
    shape1 = 2,
    shape2 = 3,
    scale = 1
  )

  out <- discretize_severity(
    x,
    step = 0.5,
    to = 2,
    method = "unbiased"
  )

  expect_identical(
    out$metadata$discretization$limited_mean_source,
    "numerical_cdf_integration"
  )

  expect_equal(
    sum(
      out$params$probs
    ),
    1,
    tolerance = 1e-10
  )
})


test_that("discretized severity remains compatible with distribution_moments", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  out <- discretize_severity(
    x,
    step = 0.25,
    to = 5,
    method = "unbiased"
  )

  moment <- distribution_moments(
    out,
    order = 1,
    type = "raw"
  )

  expect_equal(
    moment$value,
    1 - exp(-5),
    tolerance = 1e-10
  )

  expect_identical(
    moment$status,
    "finite"
  )
})


test_that("continuous custom severity can be discretized when operations are explicit", {
  custom_pdf <- function(x, params, log = FALSE) {
    out <- stats::dexp(
      x,
      rate = params$rate,
      log = log
    )

    out
  }

  custom_cdf <- function(
    q,
    params,
    lower.tail = TRUE,
    log.p = FALSE
  ) {
    stats::pexp(
      q,
      rate = params$rate,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  x <- severity_spec(
    "custom",
    rate = 1,
    custom = list(
      pdf = custom_pdf,
      cdf = custom_cdf,
      support = list(
        type = "continuous",
        lower = 0,
        upper = Inf
      )
    )
  )

  expect_error(
    discretize_severity(
      x,
      step = 1
    ),
    "`to` must be supplied"
  )

  out <- discretize_severity(
    x,
    step = 1,
    to = 4,
    method = "upper"
  )

  expect_s3_class(
    out,
    "tidyact_severity_spec"
  )

  expect_equal(
    sum(
      out$params$probs
    ),
    1,
    tolerance = 1e-12
  )
})


test_that("noncontinuous severities are not re-discretized", {
  x <- severity_spec(
    "discrete",
    values = c(
      0,
      10
    ),
    probs = c(
      0.5,
      0.5
    )
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = 10
    ),
    "continuous support"
  )

  empirical <- empirical_severity(
    c(
      10,
      20,
      30
    )
  )

  expect_error(
    discretize_severity(
      empirical,
      step = 1,
      to = 30
    ),
    "continuous support"
  )
})


test_that("discretize_severity validates x", {
  expect_error(
    discretize_severity(
      list(),
      step = 1,
      to = 10
    ),
    "`x` must be"
  )

  expect_error(
    discretize_severity(
      frequency_spec(
        "poisson",
        lambda = 2
      ),
      step = 1,
      to = 10
    ),
    "`x` must be"
  )
})


test_that("discretize_severity validates step", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  for (bad in list(
    0,
    -1,
    Inf,
    NA_real_,
    c(1, 2),
    "1"
  )) {
    expect_error(
      discretize_severity(
        x,
        step = bad,
        to = 10
      ),
      "`step`"
    )
  }
})


test_that("discretize_severity validates to and its lattice alignment", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = 0
    ),
    "`to`"
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = -1
    ),
    "`to`"
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = Inf
    ),
    "`to`"
  )

  expect_error(
    discretize_severity(
      x,
      step = 0.3,
      to = 1
    ),
    "integer multiple"
  )
})


test_that("discretize_severity validates method", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = 5,
      method = "unknown"
    )
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = 5,
      method = NA_character_
    ),
    "`method`"
  )
})


test_that("discretize_severity validates tail_prob", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  for (bad in list(
    0,
    1,
    -0.1,
    1.1,
    Inf,
    NA_real_,
    c(0.1, 0.2)
  )) {
    expect_error(
      discretize_severity(
        x,
        step = 1,
        to = 5,
        tail_prob = bad
      ),
      "`tail_prob`"
    )
  }
})


test_that("discretize_severity validates label", {
  x <- severity_spec(
    "exponential",
    rate = 1
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = 5,
      label = ""
    ),
    "`label`"
  )

  expect_error(
    discretize_severity(
      x,
      step = 1,
      to = 5,
      label = NA_character_
    ),
    "`label`"
  )
})
