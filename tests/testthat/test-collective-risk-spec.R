# Tests for collective_risk_spec --------------------------------------------

test_that("collective_risk_spec builds the expected S3 object", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 4
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 0.001
  )

  x <- collective_risk_spec(
    frequency,
    severity
  )

  expect_s3_class(
    x,
    "tidyact_collective_risk_spec"
  )

  expect_s3_class(
    x,
    "tidyact_risk_spec"
  )

  expect_identical(
    x$schema_version,
    1L
  )

  expect_identical(
    x$frequency,
    frequency
  )

  expect_identical(
    x$severity,
    severity
  )

  expect_identical(
    x$count_basis,
    "loss"
  )

  expect_null(x$label)

  expect_equal(
    x$metadata,
    list()
  )

  expect_named(
    x,
    c(
      "schema_version",
      "frequency",
      "severity",
      "count_basis",
      "label",
      "metadata"
    ),
    ignore.order = FALSE
  )
})


test_that("collective_risk_spec accepts explicit loss and payment bases", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 2
  )

  severity <- severity_spec(
    "exponential",
    rate = 0.5
  )

  loss <- collective_risk_spec(
    frequency,
    severity,
    count_basis = " LOSS "
  )

  payment <- collective_risk_spec(
    frequency,
    severity,
    count_basis = "PAYMENT"
  )

  expect_identical(
    loss$count_basis,
    "loss"
  )

  expect_identical(
    payment$count_basis,
    "payment"
  )
})


test_that("collective_risk_spec carries an optional label", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  x <- collective_risk_spec(
    frequency,
    severity,
    label = "Motor portfolio"
  )

  expect_identical(
    x$label,
    "Motor portfolio"
  )
})


test_that("count basis is inferred from component metadata", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  frequency$metadata$count_basis <- "payment"

  x <- collective_risk_spec(
    frequency,
    severity
  )

  expect_identical(
    x$count_basis,
    "payment"
  )
})


test_that("matching count-basis metadata is accepted", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  frequency$metadata$count_basis <- "payment"
  severity$metadata$count_basis <- "payment"

  x <- collective_risk_spec(
    frequency,
    severity
  )

  expect_identical(
    x$count_basis,
    "payment"
  )
})


test_that("conflicting inferred count bases are rejected", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  frequency$metadata$count_basis <- "loss"
  severity$metadata$count_basis <- "payment"

  expect_error(
    collective_risk_spec(
      frequency,
      severity
    ),
    "conflicting count-basis metadata"
  )
})


test_that("an explicit count basis resolves component metadata deliberately", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  frequency$metadata$count_basis <- "loss"
  severity$metadata$count_basis <- "payment"

  x <- collective_risk_spec(
    frequency,
    severity,
    count_basis = "loss"
  )

  expect_identical(
    x$count_basis,
    "loss"
  )
})


test_that("collective_risk_spec validates the frequency component", {
  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  expect_error(
    collective_risk_spec(
      list(),
      severity
    ),
    "`frequency` must be"
  )

  fake_frequency <- structure(
    list(
      schema_version = 1L,
      dist = "fake",
      params = list(),
      support = list(
        type = "continuous",
        lower = 0,
        upper = Inf
      ),
      data = NULL,
      custom = NULL,
      label = NULL,
      metadata = list()
    ),
    class = c(
      "tidyact_frequency_spec",
      "tidyact_distribution_spec"
    )
  )

  expect_error(
    collective_risk_spec(
      fake_frequency,
      severity
    ),
    "count support"
  )
})


test_that("collective_risk_spec validates the severity component", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  expect_error(
    collective_risk_spec(
      frequency,
      list()
    ),
    "`severity` must be"
  )

  fake_severity <- structure(
    list(
      schema_version = 1L,
      dist = "fake",
      params = list(),
      support = list(
        type = "continuous",
        lower = -1,
        upper = Inf
      ),
      data = NULL,
      custom = NULL,
      label = NULL,
      metadata = list()
    ),
    class = c(
      "tidyact_severity_spec",
      "tidyact_distribution_spec"
    )
  )

  expect_error(
    collective_risk_spec(
      frequency,
      fake_severity
    ),
    "nonnegative support"
  )
})


test_that("collective_risk_spec validates count_basis", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  expect_error(
    collective_risk_spec(
      frequency,
      severity,
      count_basis = "claims"
    ),
    "`count_basis`"
  )

  expect_error(
    collective_risk_spec(
      frequency,
      severity,
      count_basis = c(
        "loss",
        "payment"
      )
    ),
    "`count_basis`"
  )

  expect_error(
    collective_risk_spec(
      frequency,
      severity,
      count_basis = NA_character_
    ),
    "`count_basis`"
  )
})


test_that("collective_risk_spec validates label", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 3
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  expect_error(
    collective_risk_spec(
      frequency,
      severity,
      label = ""
    ),
    "`label`"
  )

  expect_error(
    collective_risk_spec(
      frequency,
      severity,
      label = "   "
    ),
    "`label`"
  )

  expect_error(
    collective_risk_spec(
      frequency,
      severity,
      label = NA_character_
    ),
    "`label`"
  )
})


test_that("empirical components can define a collective model", {
  frequency <- empirical_frequency(
    c(
      0,
      0,
      1,
      2
    )
  )

  severity <- empirical_severity(
    c(
      100,
      200,
      500
    )
  )

  x <- collective_risk_spec(
    frequency,
    severity
  )

  expect_s3_class(
    x,
    "tidyact_collective_risk_spec"
  )

  expect_identical(
    x$frequency$dist,
    "empirical"
  )

  expect_identical(
    x$severity$dist,
    "empirical"
  )
})


test_that("print method reports the collective model", {
  frequency <- frequency_spec(
    "poisson",
    lambda = 4
  )

  severity <- severity_spec(
    "gamma",
    shape = 2,
    rate = 1
  )

  x <- collective_risk_spec(
    frequency,
    severity,
    label = "Illustrative model"
  )

  expect_output(
    print(x),
    "<tidyact_collective_risk_spec>",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Frequency: Poisson",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Severity: Gamma",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Count basis: loss",
    fixed = TRUE
  )

  expect_output(
    print(x),
    "Label: Illustrative model",
    fixed = TRUE
  )

  expect_identical(
    invisible(print(x)),
    invisible(x)
  )
})
