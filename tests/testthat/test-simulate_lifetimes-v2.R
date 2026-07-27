test_that("simulate_lifetimes returns one complete sample per life", {
  lt <- tibble::tibble(
    x = 60:62,
    qx = c(0.2, 0.5, 1)
  )

  result <- simulate_lifetimes(
    data = lt,
    x = c(60, 61),
    n_sim = 50,
    seed = 123
  )

  expect_equal(nrow(result), 100L)
  expect_equal(sort(unique(result$life_id)), c(1L, 2L))

  counts <- result |>
    dplyr::count(.data$life_id)

  expect_equal(counts$n, c(50L, 50L))
  expect_equal(result$sim, result$sim_id)
  expect_equal(result$life, result$life_id)
  expect_equal(result$K, result$Kx)
  expect_equal(result$T, result$Tx)

  expect_true(
    all(
      result$Tx >= result$Kx &
        result$Tx < result$Kx + 1
    )
  )
})


test_that("simulate_lifetimes uses one stream rather than resetting per life", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  result <- simulate_lifetimes(
    data = lt,
    x = c(60, 60),
    n_sim = 100,
    seed = 123
  )

  life_1 <- result |>
    dplyr::filter(.data$life_id == 1L)

  life_2 <- result |>
    dplyr::filter(.data$life_id == 2L)

  expect_false(
    identical(life_1$Kx, life_2$Kx)
  )

  repeated <- simulate_lifetimes(
    data = lt,
    x = c(60, 60),
    n_sim = 100,
    seed = 123
  )

  expect_identical(result, repeated)
})


test_that("qx, px, and lx mortality bases are equivalent", {
  qx_table <- tibble::tibble(
    x = 60:62,
    qx = c(0.2, 0.5, 1)
  )

  px_table <- tibble::tibble(
    x = 60:62,
    px = c(0.8, 0.5, 0)
  )

  lx_table <- tibble::tibble(
    x = 60:63,
    lx = c(100, 80, 40, 0)
  )

  from_qx <- simulate_lifetimes(
    qx_table,
    x = c(60, 61),
    n_sim = 50,
    seed = 321
  )

  from_px <- simulate_lifetimes(
    px_table,
    x = c(60, 61),
    n_sim = 50,
    seed = 321
  )

  from_lx <- simulate_lifetimes(
    lx_table,
    x = c(60, 61),
    n_sim = 50,
    seed = 321
  )

  expect_equal(
    from_qx[c("life_id", "Kx", "Tx")],
    from_px[c("life_id", "Kx", "Tx")],
    tolerance = 1e-12
  )

  expect_equal(
    from_qx[c("life_id", "Kx", "Tx")],
    from_lx[c("life_id", "Kx", "Tx")],
    tolerance = 1e-12
  )
})


test_that("simulate_lifetimes canonicalizes fractional-age aliases", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  for (alias in c(
    "constant",
    "cfm",
    "constant_force"
  )) {
    result <- simulate_lifetimes(
      lt,
      x = c(60, 60),
      n_sim = 10,
      frac = alias,
      seed = 123
    )

    expect_identical(
      unique(result$frac),
      "cml"
    )
  }
})


test_that("simulate_lifetimes rejects invalid mortality instead of repairing it", {
  invalid_qx <- tibble::tibble(
    x = 60:61,
    qx = c(-0.1, 1)
  )

  invalid_px <- tibble::tibble(
    x = 60:61,
    px = c(1.1, 0)
  )

  increasing_lx <- tibble::tibble(
    x = 60:62,
    lx = c(100, 110, 0)
  )

  expect_error(
    simulate_lifetimes(
      invalid_qx,
      x = c(60, 60),
      n_sim = 10
    ),
    "between 0 and 1"
  )

  expect_error(
    simulate_lifetimes(
      invalid_px,
      x = c(60, 60),
      n_sim = 10
    ),
    "between 0 and 1"
  )

  expect_error(
    simulate_lifetimes(
      increasing_lx,
      x = c(60, 60),
      n_sim = 10
    ),
    "nonincreasing"
  )
})


test_that("simulate_lifetimes treats truncation explicitly", {
  truncated <- tibble::tibble(
    x = 60:61,
    qx = c(0.2, 0.3)
  )

  expect_warning(
    conditional <- simulate_lifetimes(
      truncated,
      x = c(60, 60),
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
    simulate_lifetimes(
      truncated,
      x = c(60, 60),
      n_sim = 20,
      truncation = "error",
      seed = 123
    ),
    "does not exhaust"
  )
})


test_that("simulate_lifetimes restores the caller random state", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  set.seed(999)
  state_before <- .Random.seed

  invisible(
    simulate_lifetimes(
      lt,
      x = c(60, 60),
      n_sim = 20,
      seed = 0
    )
  )

  expect_identical(
    .Random.seed,
    state_before
  )
})


test_that("simulate_lifetimes validates ages, simulation count, and seed", {
  lt <- tibble::tibble(
    x = 60:61,
    qx = c(0.5, 1)
  )

  expect_error(
    simulate_lifetimes(
      lt,
      x = c(60, 60.5),
      n_sim = 10
    ),
    "integer ages"
  )

  expect_error(
    simulate_lifetimes(
      lt,
      x = c(60, 60),
      n_sim = 10.5
    ),
    "positive integer"
  )

  expect_error(
    simulate_lifetimes(
      lt,
      x = c(60, 60),
      n_sim = 10,
      seed = 1.5
    ),
    "nonnegative integer"
  )
})


test_that("simulate_lifetimes requires consecutive ages through the horizon", {
  gap <- tibble::tibble(
    x = c(60, 62, 63),
    qx = c(0.2, 0.5, 1)
  )

  expect_error(
    simulate_lifetimes(
      gap,
      x = c(60, 62),
      n_sim = 10
    ),
    "consecutive"
  )
})
