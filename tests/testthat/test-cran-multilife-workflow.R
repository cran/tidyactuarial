test_that("multiple-life simulation workflow works with mc_annuity", {
  lt <- tibble::tibble(
    x = 40:100,
    qx = seq(0.002, 1, length.out = 61)
  )

  expect_no_error({
    lt |>
      simulate_lifetimes(
        x = c(60, 58),
        n_sim = 25,
        frac = "udd",
        seed = 123
      ) |>
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
  })
})
