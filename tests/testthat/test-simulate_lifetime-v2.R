test_that("simulate_lifetime reproduces a simple curtate distribution", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  result <- simulate_lifetime(
    lt,
    x = 60,
    n_sim = 10000,
    frac = "none",
    seed = 123
  )

  expect_equal(
    mean(result$Kx == 0L),
    0.5,
    tolerance = 0.02
  )

  expect_equal(
    mean(result$Kx == 1L),
    0.5,
    tolerance = 0.02
  )

  expect_true(all(is.na(result$Tx)))
  expect_false(any(result$distribution_conditioned))
})


test_that("simulate_lifetime requires x and consecutive attained ages", {
  missing_x <- tibble::tibble(
    x = 61:63,
    qx = c(0.2, 0.3, 1)
  )

  gap <- tibble::tibble(
    x = c(60, 62, 63),
    qx = c(0.2, 0.3, 1)
  )

  duplicated <- tibble::tibble(
    x = c(60, 60, 61),
    qx = c(0.2, 0.3, 1)
  )

  expect_error(
    simulate_lifetime(
      missing_x,
      x = 60,
      n_sim = 10
    ),
    "appear explicitly"
  )

  expect_error(
    simulate_lifetime(
      gap,
      x = 60,
      n_sim = 10
    ),
    "consecutive"
  )

  expect_error(
    simulate_lifetime(
      duplicated,
      x = 60,
      n_sim = 10
    ),
    "duplicated"
  )
})


test_that("simulate_lifetime supports all fractional-age assumptions", {
  lt <- tibble::tibble(
    x = 60:62,
    qx = c(0.2, 0.4, 1)
  )

  for (assumption in c("udd", "cml", "balducci")) {
    result <- simulate_lifetime(
      lt,
      x = 60,
      n_sim = 500,
      frac = assumption,
      seed = 123
    )

    expect_true(
      all(
        result$Tx >= result$Kx &
          result$Tx < result$Kx + 1
      )
    )
  }

  legacy <- simulate_lifetime(
    lt,
    x = 60,
    n_sim = 20,
    frac = "constant_force",
    seed = 123
  )

  expect_identical(
    unique(legacy$frac),
    "cml"
  )
})


test_that("simulate_lifetime handles truncated tables explicitly", {
  truncated <- tibble::tibble(
    x = 60:61,
    qx = c(0.2, 0.3)
  )

  expect_warning(
    conditional <- simulate_lifetime(
      truncated,
      x = 60,
      n_sim = 20,
      truncation = "conditional",
      seed = 123
    ),
    "conditional on death"
  )

  expect_true(
    all(conditional$distribution_conditioned)
  )

  expect_error(
    simulate_lifetime(
      truncated,
      x = 60,
      n_sim = 20,
      truncation = "error",
      seed = 123
    ),
    "does not exhaust"
  )
})


test_that("simulate_lifetime antithetic draws are paired", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  result <- simulate_lifetime(
    lt,
    x = 60,
    n_sim = 100,
    method = "antithetic",
    frac = "none",
    seed = 123
  )

  pair_matrix <- matrix(
    result$Kx,
    ncol = 2,
    byrow = TRUE
  )

  expect_true(
    all(rowSums(pair_matrix) == 1L)
  )
})


test_that("simulate_lifetime restores the caller random state", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  set.seed(999)
  state_before <- .Random.seed

  invisible(
    simulate_lifetime(
      lt,
      x = 60,
      n_sim = 20,
      seed = 0
    )
  )

  expect_identical(
    .Random.seed,
    state_before
  )
})


test_that("simulate_lifetime returns standardized aliases and distribution", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  result <- simulate_lifetime(
    lt,
    x = 60,
    n_sim = 5,
    include_distribution = TRUE,
    seed = 123
  )

  expect_equal(
    result$simulation_id,
    result$sim_id
  )

  expect_equal(
    result$curtate_lifetime,
    result$Kx
  )

  expect_equal(
    result$complete_lifetime,
    result$Tx
  )

  expect_true(
    all(
      vapply(
        result$distribution,
        inherits,
        logical(1),
        "tbl_df"
      )
    )
  )

  distribution <- result$distribution[[1]]

  expect_true(
    all(
      c(
        "prob_unconditional",
        "prob",
        "cdf"
      ) %in% names(distribution)
    )
  )

  expect_equal(
    sum(distribution$prob),
    1,
    tolerance = 1e-12
  )
})


test_that("simulate_lifetime validates age, columns, and probabilities", {
  valid <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  invalid_qx <- tibble::tibble(
    x = 60:61,
    qx = c(-0.1, 1)
  )

  expect_error(
    simulate_lifetime(
      valid,
      x = 60.5,
      n_sim = 10
    ),
    "integer age"
  )

  expect_error(
    simulate_lifetime(
      valid,
      x = 60,
      n_sim = 10,
      x_col = ""
    ),
    "nonempty"
  )

  expect_error(
    simulate_lifetime(
      invalid_qx,
      x = 60,
      n_sim = 10
    ),
    "between 0 and 1"
  )

  expect_error(
    simulate_lifetime(
      valid,
      x = 60,
      n_sim = 10,
      seed = 1.5
    ),
    "nonnegative integer"
  )
})
