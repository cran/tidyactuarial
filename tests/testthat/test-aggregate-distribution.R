# Tests for aggregate_distribution ------------------------------------------

.agg_probability_at <- function(x, loss) {
  vapply(
    loss,
    function(z) {
      idx <- which(
        abs(
          x$grid$loss -
            z
        ) <=
          1e-10 *
            max(
              1,
              abs(z)
            )
      )

      if (!length(idx)) {
        return(0)
      }

      sum(
        x$grid$probability[idx]
      )
    },
    numeric(1)
  )
}


.agg_compare_on_grid <- function(
  x,
  y,
  grid,
  tolerance = 1e-10
) {
  expect_equal(
    .agg_probability_at(
      x,
      grid
    ),
    .agg_probability_at(
      y,
      grid
    ),
    tolerance = tolerance
  )
}


.make_finite_aggregate_model <- function() {
  frequency <- frequency_spec(
    "binomial",
    size = 4,
    prob = 0.3
  )

  severity <- severity_spec(
    "discrete",
    values = c(
      0,
      10,
      20
    ),
    probs = c(
      0.2,
      0.5,
      0.3
    )
  )

  collective_risk_spec(
    frequency,
    severity
  )
}


test_that("aggregate_distribution returns the expected S3 structure", {
  model <- .make_finite_aggregate_model()

  out <- aggregate_distribution(
    model,
    method = "convolution"
  )

  expect_s3_class(
    out,
    "tidyact_aggregate_distribution"
  )

  expect_s3_class(
    out,
    "tidyact_risk_distribution"
  )

  expect_identical(
    out$schema_version,
    1L
  )

  expect_identical(
    out$method,
    "convolution"
  )

  expect_s3_class(
    out$grid,
    "tbl_df"
  )

  expect_named(
    out$grid,
    c(
      "loss",
      "probability",
      "cdf"
    ),
    ignore.order = FALSE
  )

  expect_identical(
    out$model,
    model
  )

  expect_true(
    is.list(
      out$metadata
    )
  )
})


test_that("returned probabilities and explicit tail conserve total probability", {
  model <- .make_finite_aggregate_model()

  methods <- c(
    "convolution",
    "panjer",
    "fft"
  )

  for (method in methods) {
    out <- aggregate_distribution(
      model,
      method = method
    )

    expect_true(
      all(
        out$grid$probability >= 0
      )
    )

    expect_equal(
      sum(
        out$grid$probability
      ) +
        out$tail_probability,
      1,
      tolerance = 1e-9
    )

    expect_equal(
      out$grid$cdf,
      cumsum(
        out$grid$probability
      ),
      tolerance = 1e-13
    )
  }
})


test_that("finite convolution recovers a simple exact random sum", {
  frequency <- empirical_frequency(
    2
  )

  severity <- severity_spec(
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

  model <- collective_risk_spec(
    frequency,
    severity
  )

  out <- aggregate_distribution(
    model,
    method = "convolution"
  )

  expect_equal(
    .agg_probability_at(
      out,
      c(
        0,
        10,
        20
      )
    ),
    c(
      0.25,
      0.50,
      0.25
    ),
    tolerance = 1e-13
  )

  expect_equal(
    out$tail_probability,
    0,
    tolerance = 1e-13
  )

  expect_true(
    isTRUE(
      out$metadata$exact
    )
  )
})


test_that("convolution, Panjer, and FFT agree for an exact finite model", {
  model <- .make_finite_aggregate_model()

  convolution <- aggregate_distribution(
    model,
    method = "convolution"
  )

  panjer <- aggregate_distribution(
    model,
    method = "panjer"
  )

  fft <- aggregate_distribution(
    model,
    method = "fft"
  )

  grid <- seq(
    0,
    80,
    by = 10
  )

  .agg_compare_on_grid(
    convolution,
    panjer,
    grid,
    tolerance = 1e-10
  )

  .agg_compare_on_grid(
    convolution,
    fft,
    grid,
    tolerance = 1e-9
  )

  expect_equal(
    panjer$tail_probability,
    0,
    tolerance = 1e-9
  )

  expect_equal(
    fft$tail_probability,
    0,
    tolerance = 1e-9
  )
})


test_that("exact aggregate mean agrees with aggregate_moments", {
  model <- .make_finite_aggregate_model()

  exact <- aggregate_distribution(
    model,
    method = "convolution"
  )

  grid_mean <- sum(
    exact$grid$loss *
      exact$grid$probability
  )

  actuarial_mean <- aggregate_moments(
    model,
    order = 1,
    type = "raw"
  )$value

  expect_equal(
    grid_mean,
    actuarial_mean,
    tolerance = 1e-10
  )
})


test_that("Panjer zero probability equals the frequency PGF at f0", {
  lambda <- 2

  frequency <- frequency_spec(
    "poisson",
    lambda = lambda
  )

  severity <- severity_spec(
    "discrete",
    values = c(
      0,
      10
    ),
    probs = c(
      0.2,
      0.8
    )
  )

  model <- collective_risk_spec(
    frequency,
    severity
  )

  out <- aggregate_distribution(
    model,
    method = "panjer",
    tol = 1e-10
  )

  expected_p0 <- exp(
    lambda *
      (
        0.2 -
          1
      )
  )

  expect_equal(
    .agg_probability_at(
      out,
      0
    ),
    expected_p0,
    tolerance = 1e-12
  )

  expect_lte(
    out$tail_probability,
    1e-10 +
      1e-12
  )
})


test_that("zero frequency gives an exact degenerate aggregate", {
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

  out <- aggregate_distribution(
    model,
    method = "auto"
  )

  expect_identical(
    out$method,
    "analytic"
  )

  expect_equal(
    out$grid$loss,
    0
  )

  expect_equal(
    out$grid$probability,
    1
  )

  expect_equal(
    out$tail_probability,
    0
  )

  expect_identical(
    out$metadata$analytic_case,
    "zero_frequency"
  )
})


test_that("zero severity gives an exact degenerate aggregate", {
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

  out <- aggregate_distribution(
    model,
    method = "auto"
  )

  expect_identical(
    out$method,
    "analytic"
  )

  expect_equal(
    out$grid$loss,
    0
  )

  expect_equal(
    out$grid$probability,
    1
  )

  expect_equal(
    out$tail_probability,
    0
  )
})


test_that("deterministic severity transforms count values exactly", {
  frequency <- frequency_spec(
    "binomial",
    size = 3,
    prob = 0.4
  )

  model <- collective_risk_spec(
    frequency,
    severity_spec(
      "degenerate",
      value = 5
    )
  )

  out <- aggregate_distribution(
    model,
    method = "analytic"
  )

  expected_losses <- c(
    0,
    5,
    10,
    15
  )

  expected_probs <- stats::dbinom(
    0:3,
    size = 3,
    prob = 0.4
  )

  expect_equal(
    out$grid$loss,
    expected_losses
  )

  expect_equal(
    out$grid$probability,
    expected_probs,
    tolerance = 1e-13
  )

  expect_equal(
    out$tail_probability,
    0
  )

  expect_true(
    isTRUE(
      out$metadata$exact
    )
  )
})


test_that("deterministic infinite-count model reports the truncated tail explicitly", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 4
    ),
    severity_spec(
      "degenerate",
      value = 10
    )
  )

  out <- aggregate_distribution(
    model,
    method = "analytic",
    tol = 1e-8
  )

  expect_equal(
    sum(
      out$grid$probability
    ) +
      out$tail_probability,
    1,
    tolerance = 1e-12
  )

  expect_lte(
    out$tail_probability,
    1e-8 +
      1e-12
  )

  expect_gt(
    out$tail_probability,
    0
  )
})


test_that("max_loss retains omitted probability as tail instead of renormalizing", {
  model <- .make_finite_aggregate_model()

  full <- aggregate_distribution(
    model,
    method = "convolution"
  )

  cutoff <- 30

  truncated <- aggregate_distribution(
    model,
    method = "convolution",
    max_loss = cutoff
  )

  omitted <- sum(
    full$grid$probability[
      full$grid$loss >
        cutoff
    ]
  )

  expect_equal(
    truncated$tail_probability,
    omitted,
    tolerance = 1e-12
  )

  expect_equal(
    sum(
      truncated$grid$probability
    ) +
      truncated$tail_probability,
    1,
    tolerance = 1e-12
  )

  expect_lte(
    max(
      truncated$grid$loss
    ),
    cutoff
  )
})


test_that("Panjer max_loss preserves explicit tail probability", {
  model <- .make_finite_aggregate_model()

  out <- aggregate_distribution(
    model,
    method = "panjer",
    max_loss = 30
  )

  expect_equal(
    sum(
      out$grid$probability
    ) +
      out$tail_probability,
    1,
    tolerance = 1e-10
  )

  expect_gt(
    out$tail_probability,
    0
  )

  expect_lte(
    max(
      out$grid$loss
    ),
    30
  )
})


test_that("FFT max_loss preserves explicit tail probability", {
  model <- .make_finite_aggregate_model()

  out <- aggregate_distribution(
    model,
    method = "fft",
    max_loss = 30
  )

  expect_equal(
    sum(
      out$grid$probability
    ) +
      out$tail_probability,
    1,
    tolerance = 1e-9
  )

  expect_gt(
    out$tail_probability,
    0
  )

  expect_lte(
    max(
      out$grid$loss
    ),
    30
  )
})


test_that("auto selects analytic for deterministic severity", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    severity_spec(
      "degenerate",
      value = 5
    )
  )

  out <- aggregate_distribution(
    model,
    method = "auto"
  )

  expect_identical(
    out$method,
    "analytic"
  )

  expect_identical(
    out$metadata$requested_method,
    "auto"
  )

  expect_identical(
    out$metadata$selected_method,
    "analytic"
  )
})


test_that("auto selects Panjer for a compatible lattice model", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    severity_spec(
      "discrete",
      values = c(
        0,
        10
      ),
      probs = c(
        0.2,
        0.8
      )
    )
  )

  out <- aggregate_distribution(
    model,
    method = "auto",
    tol = 1e-8
  )

  expect_identical(
    out$method,
    "panjer"
  )

  expect_identical(
    out$metadata$selected_method,
    "panjer"
  )
})


test_that("auto selects exact convolution for finite empirical frequency without Panjer", {
  frequency <- empirical_frequency(
    c(
      0,
      1,
      2,
      2
    )
  )

  severity <- severity_spec(
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

  model <- collective_risk_spec(
    frequency,
    severity
  )

  out <- aggregate_distribution(
    model,
    method = "auto"
  )

  expect_identical(
    out$method,
    "convolution"
  )

  expect_true(
    isTRUE(
      out$metadata$exact
    )
  )
})


test_that("auto can select FFT for an infinite non-Panjer custom frequency", {
  custom_pmf <- function(
    x,
    params,
    log = FALSE
  ) {
    stats::dpois(
      x,
      lambda = params$lambda,
      log = log
    )
  }

  custom_cdf <- function(
    q,
    params,
    lower.tail = TRUE,
    log.p = FALSE
  ) {
    stats::ppois(
      q,
      lambda = params$lambda,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  custom_quantile <- function(
    p,
    params,
    lower.tail = TRUE,
    log.p = FALSE
  ) {
    stats::qpois(
      p,
      lambda = params$lambda,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  frequency <- frequency_spec(
    "custom",
    lambda = 2,
    custom = list(
      pmf = custom_pmf,
      cdf = custom_cdf,
      quantile = custom_quantile,
      support = list(
        type = "count",
        lower = 0,
        upper = Inf,
        step = 1
      )
    )
  )

  severity <- severity_spec(
    "discrete",
    values = c(
      0,
      10
    ),
    probs = c(
      0.3,
      0.7
    )
  )

  model <- collective_risk_spec(
    frequency,
    severity
  )

  out <- aggregate_distribution(
    model,
    method = "auto",
    tol = 1e-8
  )

  expect_identical(
    out$method,
    "fft"
  )

  expect_lte(
    out$tail_probability,
    1e-8 +
      1e-12
  )
})


test_that("continuous severity is never discretized implicitly", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    severity_spec(
      "gamma",
      shape = 2,
      rate = 1
    )
  )

  expect_error(
    aggregate_distribution(
      model,
      method = "auto"
    ),
    "discretize_severity"
  )

  expect_error(
    aggregate_distribution(
      model,
      method = "panjer"
    ),
    "discretize_severity"
  )

  expect_error(
    aggregate_distribution(
      model,
      method = "fft"
    ),
    "discretize_severity"
  )

  expect_error(
    aggregate_distribution(
      model,
      method = "convolution"
    ),
    "finite discrete severity"
  )
})


test_that("explicit discretization enables Panjer and FFT", {
  continuous <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  lattice <- discretize_severity(
    continuous,
    step = 0.5,
    to = 8,
    method = "unbiased"
  )

  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    lattice
  )

  panjer <- aggregate_distribution(
    model,
    method = "panjer",
    tol = 1e-8
  )

  fft <- aggregate_distribution(
    model,
    method = "fft",
    tol = 1e-8
  )

  expect_identical(
    panjer$method,
    "panjer"
  )

  expect_identical(
    fft$method,
    "fft"
  )

  expect_equal(
    sum(
      panjer$grid$probability
    ) +
      panjer$tail_probability,
    1,
    tolerance = 1e-9
  )

  expect_equal(
    sum(
      fft$grid$probability
    ) +
      fft$tail_probability,
    1,
    tolerance = 1e-8
  )
})


test_that("FFT count truncation is explicit in metadata", {
  severity <- severity_spec(
    "discrete",
    values = c(
      0,
      1
    ),
    probs = c(
      0.5,
      0.5
    )
  )

  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 3
    ),
    severity
  )

  out <- aggregate_distribution(
    model,
    method = "fft",
    tol = 1e-6
  )

  expect_true(
    is.numeric(
      out$metadata$count_truncation
    )
  )

  expect_lte(
    out$metadata$count_tail_probability,
    1e-6 +
      1e-12
  )

  expect_equal(
    out$tail_probability,
    out$metadata$count_tail_probability,
    tolerance = 1e-9
  )
})


test_that("analytic method rejects unsupported nondeterministic severity", {
  model <- .make_finite_aggregate_model()

  expect_error(
    aggregate_distribution(
      model,
      method = "analytic"
    ),
    "deterministic severity"
  )
})


test_that("convolution method rejects infinite frequency support", {
  model <- collective_risk_spec(
    frequency_spec(
      "poisson",
      lambda = 2
    ),
    severity_spec(
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
  )

  expect_error(
    aggregate_distribution(
      model,
      method = "convolution"
    ),
    "finite claim-count support"
  )
})


test_that("Panjer rejects a frequency without a Panjer representation", {
  frequency <- empirical_frequency(
    c(
      0,
      1,
      2
    )
  )

  severity <- severity_spec(
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

  model <- collective_risk_spec(
    frequency,
    severity
  )

  expect_error(
    aggregate_distribution(
      model,
      method = "panjer"
    ),
    "Panjer"
  )
})


test_that("aggregate_distribution validates x", {
  expect_error(
    aggregate_distribution(
      list()
    ),
    "`x` must be"
  )

  expect_error(
    aggregate_distribution(
      frequency_spec(
        "poisson",
        lambda = 2
      )
    ),
    "`x` must be"
  )
})


test_that("aggregate_distribution validates tol", {
  model <- .make_finite_aggregate_model()

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
      aggregate_distribution(
        model,
        method = "convolution",
        tol = bad
      ),
      "`tol`"
    )
  }
})


test_that("aggregate_distribution validates max_loss", {
  model <- .make_finite_aggregate_model()

  for (bad in list(
    -1,
    Inf,
    NA_real_,
    c(10, 20),
    "10"
  )) {
    expect_error(
      aggregate_distribution(
        model,
        method = "convolution",
        max_loss = bad
      ),
      "`max_loss`"
    )
  }
})


test_that("aggregate_distribution rejects unsupported dots", {
  model <- .make_finite_aggregate_model()

  expect_error(
    aggregate_distribution(
      model,
      method = "convolution",
      something = 1
    ),
    "reserved"
  )
})


test_that("print method reports aggregate distribution information", {
  model <- .make_finite_aggregate_model()

  out <- aggregate_distribution(
    model,
    method = "convolution"
  )

  expect_output(
    print(out),
    "<tidyact_aggregate_distribution>",
    fixed = TRUE
  )

  expect_output(
    print(out),
    "Method: convolution",
    fixed = TRUE
  )

  expect_output(
    print(out),
    "Captured probability:",
    fixed = TRUE
  )

  expect_output(
    print(out),
    "Tail probability:",
    fixed = TRUE
  )

  expect_identical(
    invisible(
      print(out)
    ),
    invisible(out)
  )
})
