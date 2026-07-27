test_that("inverse simulation is stable when the CDF slightly exceeds one", {
  lt <- tibble::tibble(
    x = 40:100,
    qx = seq(0.002, 1, length.out = 61)
  )

  expect_no_error(
    simulated <- simulate_lifetimes(
      data = lt,
      x = c(60, 58),
      n_sim = 25,
      frac = "udd",
      seed = 123
    )
  )

  expect_no_error(
    result <- simulated |>
      mc_multilife_status(status = "joint") |>
      mc_annuity(
        i = 0.04,
        type = "whole",
        payment = 1,
        k = 1,
        timing = "due",
        col_K = "K_status",
        col_T = "T_status"
      )
  )

  expect_equal(nrow(result), 25L)
  expect_true(all(is.finite(result$pv_annuity)))
})


test_that("simulate_lifetime returns a nondecreasing bounded CDF", {
  lt <- tibble::tibble(
    x = 40:100,
    qx = seq(0.002, 1, length.out = 61)
  )

  result <- simulate_lifetime(
    lt = lt,
    x = 58,
    n_sim = 5,
    include_distribution = TRUE,
    seed = 123
  )

  cdf <- result$distribution[[1]]$cdf

  expect_false(anyNA(cdf))
  expect_true(all(diff(cdf) >= 0))
  expect_true(all(cdf >= 0 & cdf <= 1))
  expect_equal(tail(cdf, 1), 1)
})
