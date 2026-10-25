# Tests for aggregate_moments ------------------------------------------------

make_compound_poisson_gamma <- function(
  lambda = 4,
  shape = 2,
  rate = 0.5
) {
  collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = lambda
    ),
    severity_spec(
      "gamma",
      shape = shape,
      rate = rate
    )
  )
}


test_that("aggregate_moments returns the expected tibble schema", {
  model <- make_compound_poisson_gamma()

  out <- aggregate_moments(
    model,
    order = 1:4
  )

  expect_s3_class(
    out,
    "tbl_df"
  )

  expect_named(
    out,
    c(
      "type",
      "order",
      "value",
      "status"
    ),
    ignore.order = FALSE
  )

  expect_identical(
    out$type,
    rep(
      "raw",
      4
    )
  )

  expect_identical(
    out$order,
    1:4
  )

  expect_true(
    all(
      out$status == "finite"
    )
  )
})


test_that("compound Poisson raw moments satisfy analytical identities", {
  lambda <- 3
  shape <- 2
  rate <- 0.5

  model <- make_compound_poisson_gamma(
    lambda = lambda,
    shape = shape,
    rate = rate
  )

  mu1 <- shape / rate
  mu2 <- shape * (shape + 1) / rate^2
  mu3 <- shape * (shape + 1) * (shape + 2) / rate^3
  mu4 <- shape * (shape + 1) * (shape + 2) *
    (shape + 3) / rate^4

  expected <- c(
    lambda * mu1,
    lambda * mu2 +
      lambda^2 * mu1^2,
    lambda * mu3 +
      3 * lambda^2 * mu1 * mu2 +
      lambda^3 * mu1^3,
    lambda * mu4 +
      4 * lambda^2 * mu1 * mu3 +
      3 * lambda^2 * mu2^2 +
      6 * lambda^3 * mu1^2 * mu2 +
      lambda^4 * mu1^4
  )

  out <- aggregate_moments(
    model,
    order = 1:4,
    type = "raw"
  )

  expect_equal(
    out$value,
    expected,
    tolerance = 1e-10
  )

  expect_true(
    all(
      out$status == "finite"
    )
  )
})


test_that("compound Poisson cumulants equal lambda times severity raw moments", {
  lambda <- 5
  shape <- 3
  rate <- 2

  model <- make_compound_poisson_gamma(
    lambda = lambda,
    shape = shape,
    rate = rate
  )

  severity_raw <- c(
    shape / rate,
    shape * (shape + 1) / rate^2,
    shape * (shape + 1) * (shape + 2) / rate^3,
    shape * (shape + 1) * (shape + 2) *
      (shape + 3) / rate^4
  )

  out <- aggregate_moments(
    model,
    order = 1:4,
    type = "cumulant"
  )

  expect_equal(
    out$value,
    lambda * severity_raw,
    tolerance = 1e-11
  )

  expect_true(
    all(
      out$status == "finite"
    )
  )
})


test_that("compound Poisson variance equals lambda E[X^2]", {
  lambda <- 4
  shape <- 2.5
  rate <- 1.2

  model <- make_compound_poisson_gamma(
    lambda = lambda,
    shape = shape,
    rate = rate
  )

  raw2_x <- shape * (shape + 1) / rate^2

  out <- aggregate_moments(
    model,
    order = 2,
    type = "central"
  )

  expect_equal(
    out$value,
    lambda * raw2_x,
    tolerance = 1e-11
  )

  expect_identical(
    out$status,
    "finite"
  )
})


test_that("general random-sum mean and variance identities are recovered", {
  m <- 20
  p <- 0.3
  shape <- 2
  rate <- 0.5

  frequency <- frequency_spec(
    "binomial",
    size = m,
    prob = p
  )

  severity <- severity_spec(
    "gamma",
    shape = shape,
    rate = rate
  )

  model <- collective_risk_spec(
    frequency,
    severity
  )

  mean_n <- m * p
  var_n <- m * p * (1 - p)

  mean_x <- shape / rate
  var_x <- shape / rate^2

  expected_mean <- mean_n * mean_x

  expected_var <- mean_n * var_x +
    var_n * mean_x^2

  raw <- aggregate_moments(
    model,
    order = 1,
    type = "raw"
  )

  central <- aggregate_moments(
    model,
    order = 2,
    type = "central"
  )

  expect_equal(
    raw$value,
    expected_mean,
    tolerance = 1e-11
  )

  expect_equal(
    central$value,
    expected_var,
    tolerance = 1e-10
  )
})


test_that("deterministic severity scales frequency moments exactly", {
  c_value <- 7

  frequency <- frequency_spec(
    "binomial",
    size = 12,
    prob = 0.35
  )

  severity <- severity_spec(
    "degenerate",
    value = c_value
  )

  model <- collective_risk_spec(
    frequency,
    severity
  )

  aggregate <- aggregate_moments(
    model,
    order = 1:4,
    type = "raw"
  )

  frequency_raw <- distribution_moments(
    frequency,
    order = 1:4,
    type = "raw"
  )

  expect_equal(
    aggregate$value,
    c_value^(1:4) *
      frequency_raw$value,
    tolerance = 1e-10
  )
})


test_that("generic cumulants work for a scaled binomial frequency", {
  m <- 15
  p <- 0.2
  q <- 1 - p
  c_value <- 4

  model <- collective_risk_spec(
    frequency_spec(
      "binomial",
      size = m,
      prob = p
    ),
    severity_spec(
      "degenerate",
      value = c_value
    )
  )

  out <- aggregate_moments(
    model,
    order = 1:3,
    type = "cumulant"
  )

  expected <- c(
    c_value * m * p,
    c_value^2 * m * p * q,
    c_value^3 * m * p * q * (1 - 2 * p)
  )

  expect_equal(
    out$value,
    expected,
    tolerance = 1e-10
  )

  expect_true(
    all(
      out$status == "finite"
    )
  )
})


test_that("zero claim frequency makes the aggregate degenerate at zero", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 0
    ),
    severity_spec(
      "gamma",
      shape = 2,
      rate = 1
    )
  )

  raw <- aggregate_moments(
    model,
    order = 1:6,
    type = "raw"
  )

  central <- aggregate_moments(
    model,
    order = 1:6,
    type = "central"
  )

  cumulant <- aggregate_moments(
    model,
    order = 1:6,
    type = "cumulant"
  )

  expect_equal(
    raw$value,
    rep(
      0,
      6
    )
  )

  expect_equal(
    central$value,
    rep(
      0,
      6
    )
  )

  expect_equal(
    cumulant$value,
    rep(
      0,
      6
    )
  )

  expect_true(
    all(
      raw$status == "finite"
    )
  )
})


test_that("zero severity makes the aggregate degenerate at zero", {
  model <- collective_risk_spec(
    frequency_spec(
      "negative_binomial",
      size = 2,
      prob = 0.4
    ),
    severity_spec(
      "degenerate",
      value = 0
    )
  )

  out <- aggregate_moments(
    model,
    order = 1:5,
    type = "raw"
  )

  expect_equal(
    out$value,
    rep(
      0,
      5
    )
  )

  expect_true(
    all(
      out$status == "finite"
    )
  )
})


test_that("zeroth raw and central aggregate moments equal one", {
  model <- make_compound_poisson_gamma()

  raw <- aggregate_moments(
    model,
    order = 0,
    type = "raw"
  )

  central <- aggregate_moments(
    model,
    order = 0,
    type = "central"
  )

  expect_equal(
    raw$value,
    1
  )

  expect_equal(
    central$value,
    1
  )

  expect_identical(
    raw$status,
    "finite"
  )

  expect_identical(
    central$status,
    "finite"
  )
})


test_that("first central aggregate moment equals zero", {
  model <- make_compound_poisson_gamma()

  out <- aggregate_moments(
    model,
    order = 1,
    type = "central"
  )

  expect_equal(
    out$value,
    0,
    tolerance = 1e-14
  )

  expect_identical(
    out$status,
    "finite"
  )
})


test_that("requested order and duplicates are preserved", {
  model <- make_compound_poisson_gamma()

  out <- aggregate_moments(
    model,
    order = c(
      3,
      1,
      2,
      1
    ),
    type = "raw"
  )

  expect_identical(
    out$order,
    c(
      3L,
      1L,
      2L,
      1L
    )
  )

  expect_equal(
    out$value[[2L]],
    out$value[[4L]]
  )
})


test_that("heavy-tailed severity statuses propagate without approximation", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    severity_spec(
      "pareto",
      shape = 1.5,
      scale = 1000
    )
  )

  raw <- aggregate_moments(
    model,
    order = 1:2,
    type = "raw"
  )

  expect_identical(
    raw$status,
    c(
      "finite",
      "infinite"
    )
  )

  central <- aggregate_moments(
    model,
    order = 2,
    type = "central"
  )

  expect_identical(
    central$status,
    "infinite"
  )

  cumulant <- aggregate_moments(
    model,
    order = 1:2,
    type = "cumulant"
  )

  expect_identical(
    cumulant$status,
    c(
      "finite",
      "undefined"
    )
  )
})


test_that("an infinite aggregate mean makes central moments undefined", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    severity_spec(
      "pareto",
      shape = 0.8,
      scale = 1000
    )
  )

  central <- aggregate_moments(
    model,
    order = 1:3,
    type = "central"
  )

  expect_true(
    all(
      central$status == "undefined"
    )
  )
})


test_that("empirical components integrate with aggregate_moments", {
  frequency <- empirical_frequency(
    c(
      0,
      1
    )
  )

  severity <- empirical_severity(
    c(
      10,
      20
    )
  )

  model <- collective_risk_spec(
    frequency,
    severity
  )

  # E[N] = 0.5 and E[X] = 15.
  mean_out <- aggregate_moments(
    model,
    order = 1,
    type = "raw"
  )

  expect_equal(
    mean_out$value,
    7.5,
    tolerance = 1e-12
  )

  # Var(S) = E[N] Var(X) + Var(N) E[X]^2
  # = 0.5 * 25 + 0.25 * 225 = 68.75.
  var_out <- aggregate_moments(
    model,
    order = 2,
    type = "central"
  )

  expect_equal(
    var_out$value,
    68.75,
    tolerance = 1e-12
  )
})


test_that("aggregate_moments validates x", {
  expect_error(
    aggregate_moments(
      frequency_spec(
        "poisson",
        lambda = 2
      )
    ),
    "`x` must be"
  )

  expect_error(
    aggregate_moments(
      list()
    ),
    "`x` must be"
  )
})


test_that("aggregate_moments validates moment order", {
  model <- make_compound_poisson_gamma()

  expect_error(
    aggregate_moments(
      model,
      order = numeric()
    ),
    "nonempty vector"
  )

  expect_error(
    aggregate_moments(
      model,
      order = -1
    ),
    "nonnegative integers"
  )

  expect_error(
    aggregate_moments(
      model,
      order = 1.5
    ),
    "nonnegative integers"
  )

  expect_error(
    aggregate_moments(
      model,
      order = NA_real_
    ),
    "nonnegative integers"
  )
})


test_that("aggregate_moments validates type", {
  model <- make_compound_poisson_gamma()

  expect_error(
    aggregate_moments(
      model,
      type = "factorial"
    ),
    "`type`"
  )

  expect_error(
    aggregate_moments(
      model,
      type = c(
        "raw",
        "central"
      )
    ),
    "`type`"
  )
})


test_that("zeroth cumulant is rejected", {
  model <- make_compound_poisson_gamma()

  expect_error(
    aggregate_moments(
      model,
      order = 0,
      type = "cumulant"
    ),
    "positive integers"
  )
})
