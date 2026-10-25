# Tests for distribution_convolution ----------------------------------------

test_that("pairwise convolution creates a lazy distribution specification", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  y <- severity_spec(
    "lognormal",
    meanlog = 1,
    sdlog = 0.5
  )

  z <- distribution_convolution(
    x,
    y
  )

  expect_s3_class(
    z,
    "tidyact_convolution_spec"
  )

  expect_s3_class(
    z,
    "tidyact_distribution_spec"
  )

  expect_identical(
    z$schema_version,
    1L
  )

  expect_identical(
    z$dist,
    "convolution"
  )

  expect_identical(
    z$params$mode,
    "pair"
  )

  expect_identical(
    z$params$components[[1L]],
    x
  )

  expect_identical(
    z$params$components[[2L]],
    y
  )

  expect_null(z$data)
  expect_null(z$custom)

  expect_equal(
    z$metadata,
    list()
  )
})


test_that("iid convolution stores the component and number of copies", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  z <- distribution_convolution(
    x,
    n = 4
  )

  expect_s3_class(
    z,
    "tidyact_convolution_spec"
  )

  expect_identical(
    z$params$mode,
    "iid"
  )

  expect_identical(
    z$params$component,
    x
  )

  expect_identical(
    z$params$n,
    4L
  )
})


test_that("one-fold convolution returns the original specification", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1,
    label = "Original"
  )

  z <- distribution_convolution(
    x,
    n = 1
  )

  expect_identical(
    z,
    x
  )
})


test_that("one-fold convolution may change only the label of a copy", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1,
    label = "Original"
  )

  z <- distribution_convolution(
    x,
    n = 1,
    label = "Renamed"
  )

  expect_identical(
    z$dist,
    x$dist
  )

  expect_identical(
    z$params,
    x$params
  )

  expect_identical(
    z$support,
    x$support
  )

  expect_identical(
    z$label,
    "Renamed"
  )

  expect_identical(
    x$label,
    "Original"
  )
})


test_that("zero-fold convolution is degenerate at zero", {
  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  z_severity <- distribution_convolution(
    severity,
    n = 0
  )

  expect_identical(
    z_severity$support$type,
    "discrete"
  )

  expect_equal(
    z_severity$support$lower,
    0
  )

  expect_equal(
    z_severity$support$upper,
    0
  )

  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  z_frequency <- distribution_convolution(
    frequency,
    n = 0
  )

  expect_identical(
    z_frequency$support$type,
    "count"
  )

  expect_equal(
    z_frequency$support$lower,
    0
  )

  expect_equal(
    z_frequency$support$upper,
    0
  )

  expect_equal(
    z_frequency$support$step,
    1
  )
})


test_that("continuous pairwise support bounds are propagated", {
  x <- severity_spec(
    "degenerate",
    value = 5
  )

  y <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  z <- distribution_convolution(
    x,
    y
  )

  expect_identical(
    z$support$type,
    "continuous"
  )

  expect_equal(
    z$support$lower,
    5
  )

  expect_equal(
    z$support$upper,
    Inf
  )
})


test_that("finite discrete support bounds are added", {
  x <- severity_spec(
    "discrete",
    values = c(
      0,
      10
    ),
    probs = c(
      0.4,
      0.6
    )
  )

  y <- severity_spec(
    "discrete",
    values = c(
      5,
      20
    ),
    probs = c(
      0.7,
      0.3
    )
  )

  z <- distribution_convolution(
    x,
    y
  )

  expect_identical(
    z$support$type,
    "discrete"
  )

  expect_equal(
    z$support$lower,
    5
  )

  expect_equal(
    z$support$upper,
    30
  )
})


test_that("count plus count remains count", {
  x <- frequency_spec(
    "binomial",
    size = 10,
    prob = 0.3
  )

  y <- frequency_spec(
    "binomial",
    size = 7,
    prob = 0.4
  )

  z <- distribution_convolution(
    x,
    y
  )

  expect_identical(
    z$support$type,
    "count"
  )

  expect_equal(
    z$support$lower,
    0
  )

  expect_equal(
    z$support$upper,
    17
  )

  expect_equal(
    z$support$step,
    1
  )
})


test_that("iid count convolution preserves count step and scales bounds", {
  x <- frequency_spec(
    "binomial",
    size = 7,
    prob = 0.4
  )

  z <- distribution_convolution(
    x,
    n = 5
  )

  expect_identical(
    z$support$type,
    "count"
  )

  expect_equal(
    z$support$lower,
    0
  )

  expect_equal(
    z$support$upper,
    35
  )

  expect_equal(
    z$support$step,
    1
  )
})


test_that("continuous plus count is classified as continuous", {
  continuous <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  count <- frequency_spec(
    "binomial",
    size = 4,
    prob = 0.3
  )

  z <- distribution_convolution(
    continuous,
    count
  )

  expect_identical(
    z$support$type,
    "continuous"
  )

  expect_equal(
    z$support$lower,
    0
  )

  expect_equal(
    z$support$upper,
    Inf
  )
})


test_that("iid support bounds scale algebraically", {
  x <- severity_spec(
    "discrete",
    values = c(
      2,
      5
    ),
    probs = c(
      0.4,
      0.6
    )
  )

  z <- distribution_convolution(
    x,
    n = 4
  )

  expect_identical(
    z$support$type,
    "discrete"
  )

  expect_equal(
    z$support$lower,
    8
  )

  expect_equal(
    z$support$upper,
    20
  )
})


test_that("nested convolution specifications remain composable", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  y <- severity_spec(
    "lognormal",
    meanlog = 1,
    sdlog = 0.5
  )

  first <- distribution_convolution(
    x,
    y
  )

  second <- distribution_convolution(
    first,
    x
  )

  expect_s3_class(
    second,
    "tidyact_convolution_spec"
  )

  expect_identical(
    second$params$components[[1L]],
    first
  )

  expect_identical(
    second$support$type,
    "continuous"
  )
})


test_that("convolution specification remains lazy", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  z <- distribution_convolution(
    x,
    n = 3
  )

  forbidden <- c(
    "pmf",
    "pdf",
    "cdf",
    "quantile",
    "probabilities",
    "density",
    "algorithm",
    "method",
    "fft",
    "panjer",
    "step_size",
    "tol"
  )

  expect_false(
    any(
      forbidden %in% names(z)
    )
  )

  expect_equal(
    z$metadata,
    list()
  )
})


test_that("convolution label is stored without altering components", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1,
    label = "X"
  )

  y <- severity_spec(
    "lognormal",
    meanlog = 1,
    sdlog = 0.5,
    label = "Y"
  )

  z <- distribution_convolution(
    x,
    y,
    label = "X plus Y"
  )

  expect_identical(
    z$label,
    "X plus Y"
  )

  expect_identical(
    x$label,
    "X"
  )

  expect_identical(
    y$label,
    "Y"
  )
})


test_that("distribution_convolution requires exactly one of y and n", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  y <- severity_spec(
    "lognormal",
    meanlog = 1,
    sdlog = 0.5
  )

  expect_error(
    distribution_convolution(x),
    "exactly one"
  )

  expect_error(
    distribution_convolution(
      x,
      y = y,
      n = 2
    ),
    "exactly one"
  )
})


test_that("distribution_convolution validates distribution inputs", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  expect_error(
    distribution_convolution(
      list(),
      n = 2
    ),
    "`x` must be"
  )

  expect_error(
    distribution_convolution(
      x,
      y = list()
    ),
    "`y` must be"
  )
})


test_that("distribution_convolution validates n", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  expect_error(
    distribution_convolution(
      x,
      n = -1
    ),
    "nonnegative integer"
  )

  expect_error(
    distribution_convolution(
      x,
      n = 1.5
    ),
    "nonnegative integer"
  )

  expect_error(
    distribution_convolution(
      x,
      n = NA_real_
    ),
    "nonnegative integer"
  )

  expect_error(
    distribution_convolution(
      x,
      n = Inf
    ),
    "nonnegative integer"
  )

  expect_error(
    distribution_convolution(
      x,
      n = c(
        2,
        3
      )
    ),
    "nonnegative integer"
  )
})


test_that("distribution_convolution validates label", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  expect_error(
    distribution_convolution(
      x,
      n = 2,
      label = ""
    ),
    "`label`"
  )

  expect_error(
    distribution_convolution(
      x,
      n = 2,
      label = "   "
    ),
    "`label`"
  )

  expect_error(
    distribution_convolution(
      x,
      n = 2,
      label = NA_character_
    ),
    "`label`"
  )
})


test_that("print method describes pairwise convolution", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1,
    label = "X"
  )

  y <- severity_spec(
    "lognormal",
    meanlog = 1,
    sdlog = 0.5,
    label = "Y"
  )

  z <- distribution_convolution(
    x,
    y,
    label = "Combined"
  )

  expect_output(
    print(z),
    "<tidyact_convolution_spec>",
    fixed = TRUE
  )

  expect_output(
    print(z),
    "Mode: pairwise",
    fixed = TRUE
  )

  expect_output(
    print(z),
    "\"X\" + \"Y\"",
    fixed = TRUE
  )

  expect_output(
    print(z),
    "Label: Combined",
    fixed = TRUE
  )

  expect_identical(
    invisible(print(z)),
    invisible(z)
  )
})


test_that("print method describes iid convolution", {
  x <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  z <- distribution_convolution(
    x,
    n = 3
  )

  expect_output(
    print(z),
    "Mode: iid",
    fixed = TRUE
  )

  expect_output(
    print(z),
    "Copies: 3",
    fixed = TRUE
  )
})
