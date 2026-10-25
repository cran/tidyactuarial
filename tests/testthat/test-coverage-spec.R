test_that("coverage_spec() records its default contract without extra state", {
  x <- coverage_spec()

  expect_s3_class(x, "tidyact_coverage_spec")
  expect_type(x, "list")
  expect_identical(
    names(x),
    c(
      "deductible", "limit", "coinsurance", "inflation_factor",
      "deductible_type", "index_deductible", "index_limit", "label"
    )
  )
  expect_identical(x$deductible, 0)
  expect_identical(x$limit, Inf)
  expect_identical(x$coinsurance, 1)
  expect_identical(x$inflation_factor, 1)
  expect_identical(x$deductible_type, "ordinary")
  expect_identical(x$index_deductible, FALSE)
  expect_identical(x$index_limit, FALSE)
  expect_null(x$label)
})

test_that("coverage_spec() preserves explicit scalar contract terms", {
  x <- coverage_spec(
    deductible = 500,
    limit = 10000,
    coinsurance = 0.8,
    inflation_factor = 1.1,
    deductible_type = "franchise",
    index_deductible = TRUE,
    index_limit = FALSE,
    label = "Coverage A"
  )

  expect_identical(x$deductible, 500)
  expect_identical(x$limit, 10000)
  expect_identical(x$coinsurance, 0.8)
  expect_identical(x$inflation_factor, 1.1)
  expect_identical(x$deductible_type, "franchise")
  expect_true(x$index_deductible)
  expect_false(x$index_limit)
  expect_identical(x$label, "Coverage A")
})

test_that("coverage_spec() permits mathematically degenerate coverage", {
  x <- coverage_spec(
    deductible = 0L,
    limit = 0L,
    coinsurance = 0L,
    inflation_factor = 1L,
    index_deductible = TRUE,
    index_limit = TRUE
  )

  expect_identical(x$deductible, 0)
  expect_identical(x$limit, 0)
  expect_identical(x$coinsurance, 0)
  expect_identical(x$inflation_factor, 1)
  expect_true(x$index_deductible)
  expect_true(x$index_limit)

  inf <- coverage_spec(limit = Inf, index_limit = TRUE, inflation_factor = 2)
  expect_identical(inf$limit, Inf)
  expect_true(inf$index_limit)
})

test_that("coverage_spec() does not impose a ground-up loss limit convention", {
  x <- coverage_spec(deductible = 1000, limit = 100)
  expect_identical(x$deductible, 1000)
  expect_identical(x$limit, 100)
  # A cap on the post-deductible layer need not exceed the deductible.
})

test_that("coverage_spec() rejects invalid deductibles", {
  invalid <- list(-1, Inf, -Inf, NA_real_, NaN, "10", c(1, 2), NULL, TRUE)
  for (value in invalid) {
    expect_error(coverage_spec(deductible = value), "`deductible`")
  }
})

test_that("coverage_spec() rejects invalid limits", {
  invalid <- list(-1, -Inf, NA_real_, NaN, "100", c(1, 2), NULL, TRUE)
  for (value in invalid) {
    expect_error(coverage_spec(limit = value), "`limit`")
  }
})

test_that("coverage_spec() bounds coinsurance between zero and one", {
  expect_no_error(coverage_spec(coinsurance = 0))
  expect_no_error(coverage_spec(coinsurance = 1))

  invalid <- list(-0.01, 1.01, Inf, -Inf, NA_real_, NaN,
                  "0.8", c(0.5, 0.8), NULL, TRUE)
  for (value in invalid) {
    expect_error(coverage_spec(coinsurance = value), "`coinsurance`")
  }
})

test_that("coverage_spec() requires a strictly positive finite inflation factor", {
  expect_no_error(coverage_spec(inflation_factor = 0.5))
  expect_no_error(coverage_spec(inflation_factor = 1))
  expect_no_error(coverage_spec(inflation_factor = 1.2))

  invalid <- list(0, -1, Inf, -Inf, NA_real_, NaN,
                  "1.1", c(1, 1.1), NULL, TRUE)
  for (value in invalid) {
    expect_error(coverage_spec(inflation_factor = value), "`inflation_factor`")
  }
})

test_that("coverage_spec() validates exact deductible-type names", {
  expect_identical(coverage_spec()$deductible_type, "ordinary")
  expect_identical(
    coverage_spec(deductible_type = "ordinary")$deductible_type,
    "ordinary"
  )
  expect_identical(
    coverage_spec(deductible_type = "franchise")$deductible_type,
    "franchise"
  )

  invalid <- list("ord", "Franchise", "", NA_character_,
                  1, NULL, c("ordinary", "other"))
  for (value in invalid) {
    expect_error(coverage_spec(deductible_type = value), "`deductible_type`")
  }
})

test_that("coverage_spec() requires true logical scalar indexation flags", {
  invalid <- list(1, 0, NA, NULL, "TRUE", c(TRUE, FALSE))
  for (value in invalid) {
    expect_error(coverage_spec(index_deductible = value), "`index_deductible`")
    expect_error(coverage_spec(index_limit = value), "`index_limit`")
  }
  expect_true(coverage_spec(index_deductible = TRUE)$index_deductible)
  expect_true(coverage_spec(index_limit = TRUE)$index_limit)
})

test_that("coverage_spec() validates optional contract labels", {
  expect_null(coverage_spec(label = NULL)$label)
  expect_identical(coverage_spec(label = "  Contract A ")$label,
                   "  Contract A ")
  expect_identical(coverage_spec(label = "Cobertura básica")$label,
                   "Cobertura básica")

  invalid <- list("", "   ", NA_character_, 123, TRUE,
                  character(), c("A", "B"), factor("A"))
  for (value in invalid) {
    expect_error(coverage_spec(label = value), "`label`")
  }
})

test_that("coverage_spec() only records contract terms", {
  x <- coverage_spec(deductible = 300, limit = 1000)
  expect_false(any(c(
    "frequency", "severity", "losses", "payments", "pmf", "cdf",
    "mean", "variance", "aggregate_loss", "n_sim"
  ) %in% names(x)))
  expect_identical(class(x), "tidyact_coverage_spec")
})
