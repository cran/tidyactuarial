test_that("summary_mc reproduces elementary Monte Carlo summary statistics", {
  sim <- tibble::tibble(
    present_value = 1:5
  )

  out <- summary_mc(
    sim,
    probs = c(0, 0.5, 1),
    var_probs = c(0.5, 0.8)
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)

  expect_equal(out$n_sim, 5L)
  expect_equal(out$n_total, 5L)
  expect_equal(out$n_missing, 0L)

  expect_equal(out$mean, 3)
  expect_equal(out$variance, 2.5)
  expect_equal(out$sd, sqrt(2.5))
  expect_equal(out$se_mean, sqrt(2.5) / sqrt(5))
  expect_equal(out$min, 1)
  expect_equal(out$max, 5)

  expect_equal(out$q000, 1)
  expect_equal(out$q500, 3)
  expect_equal(out$q1000, 5)

  expect_equal(out$VaR_500, 3)
  expect_equal(out$VaR_800, 4)

  # TVaR includes only the required fraction of any VaR atom.
  expect_equal(out$TVaR_500, 4.2)
  expect_equal(out$TVaR_800, 5)
})

test_that("summary_mc groups explicitly by actuarial scenario columns", {
  sim <- tibble::tibble(
    scenario = c("A", "A", "B", "B"),
    loss = c(1, 3, 2, 4)
  )

  out <- summary_mc(
    sim,
    col_value = "loss",
    by = "scenario",
    probs = 0.5,
    var_probs = 0.5
  )

  expect_equal(out$scenario, c("A", "B"))
  expect_equal(out$n_sim, c(2L, 2L))
  expect_equal(out$mean, c(2, 3))
  expect_equal(out$q500, c(2, 3))
  expect_equal(out$VaR_500, c(1, 2))
  expect_equal(out$TVaR_500, c(3, 4))
})

test_that("summary_mc respects an existing dplyr grouping structure", {
  sim <- tibble::tibble(
    t = c(0, 0, 1, 1),
    reserve = c(10, 14, 20, 24)
  ) |>
    dplyr::group_by(t)

  out <- summary_mc(
    sim,
    col_value = "reserve",
    probs = 0.5,
    var_probs = 0.5
  )

  expect_equal(out$t, c(0, 1))
  expect_equal(out$mean, c(12, 22))
})

test_that("summary_mc removes missing values only when requested", {
  sim <- tibble::tibble(
    present_value = c(10, NA, 20, 30)
  )

  removed <- summary_mc(
    sim,
    probs = 0.5,
    var_probs = 0.5,
    na_rm = TRUE
  )

  expect_equal(removed$n_total, 4L)
  expect_equal(removed$n_missing, 1L)
  expect_equal(removed$n_sim, 3L)
  expect_equal(removed$mean, 20)
  expect_equal(removed$q500, 20)

  kept <- summary_mc(
    sim,
    probs = 0.5,
    var_probs = 0.5,
    na_rm = FALSE
  )

  expect_equal(kept$n_total, 4L)
  expect_equal(kept$n_missing, 1L)
  expect_equal(kept$n_sim, 4L)
  expect_true(is.na(kept$mean))
  expect_true(is.na(kept$variance))
  expect_true(is.na(kept$q500))
  expect_true(is.na(kept$VaR_500))
  expect_true(is.na(kept$TVaR_500))
})

test_that("summary_mc handles an empty effective sample and a one-value sample", {
  empty_effective <- summary_mc(
    tibble::tibble(present_value = c(NA_real_, NA_real_)),
    probs = 0.5,
    var_probs = 0.5,
    na_rm = TRUE
  )

  expect_equal(empty_effective$n_sim, 0L)
  expect_true(is.na(empty_effective$mean))
  expect_true(is.na(empty_effective$q500))

  singleton <- summary_mc(
    tibble::tibble(present_value = 100),
    probs = 0.5,
    var_probs = 0.5
  )

  expect_equal(singleton$n_sim, 1L)
  expect_equal(singleton$mean, 100)
  expect_true(is.na(singleton$variance))
  expect_true(is.na(singleton$sd))
  expect_true(is.na(singleton$se_mean))
  expect_equal(singleton$q500, 100)
  expect_equal(singleton$VaR_500, 100)
  expect_equal(singleton$TVaR_500, 100)
})

test_that("summary_mc transitional argument names reproduce the current API", {
  sim <- tibble::tibble(
    group = c("A", "A", "B", "B"),
    value = c(2, 4, 6, 8)
  )

  current <- summary_mc(
    sim,
    col_value = "value",
    by = "group",
    probs = 0.5,
    var_probs = 0.5
  )

  legacy <- summary_mc(
    data = sim,
    value_col = "value",
    group_cols = "group",
    probs = 0.5,
    var_probs = 0.5
  )

  expect_equal(legacy, current)
})

test_that("summary_mc validates its public column interface", {
  sim <- tibble::tibble(
    value = c(1, 2, 3),
    label = c("a", "b", "c")
  )

  expect_error(
    summary_mc(sim, col_value = "missing"),
    "not found"
  )

  expect_error(
    summary_mc(sim, col_value = "label"),
    "numeric column"
  )

  expect_error(
    summary_mc(sim, col_value = "value", by = "missing"),
    "grouping columns"
  )

  expect_error(
    summary_mc(sim, col_value = "value", probs = 1.2),
    "between 0 and 1"
  )

  expect_error(
    summary_mc(sim, col_value = "value", var_probs = -0.1),
    "between 0 and 1"
  )
})
