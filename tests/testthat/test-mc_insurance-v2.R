test_that("mc_insurance whole-life end-of-year present values are correct", {
  simulations <- tibble::tibble(
    Kx = c(0, 2, 4)
  )

  result <- mc_insurance(
    simulations,
    i = 0.05,
    benefit = 100,
    type = "whole",
    timing = "end_of_year"
  )

  v <- 1 / 1.05

  expect_equal(result$benefit_time, c(1, 3, 5))
  expect_equal(
    result$pv_benefit,
    100 * v^c(1, 3, 5),
    tolerance = 1e-12
  )
})


test_that("mc_insurance discrete term and pure-endowment boundaries are correct", {
  simulations <- tibble::tibble(
    Kx = c(0, 1, 2, 3)
  )

  term_result <- mc_insurance(
    simulations,
    i = 0,
    benefit = 1000,
    type = "term",
    n = 3,
    timing = "end_of_year"
  )

  pure_result <- mc_insurance(
    simulations,
    i = 0,
    benefit = 1000,
    type = "pure_endowment",
    n = 3,
    timing = "end_of_year"
  )

  expect_equal(term_result$benefit_indicator, c(1, 1, 1, 0))
  expect_equal(term_result$pv_benefit, c(1000, 1000, 1000, 0))

  expect_equal(pure_result$benefit_indicator, c(0, 0, 0, 1))
  expect_equal(pure_result$pv_benefit, c(0, 0, 0, 1000))
})


test_that("mc_insurance endowment pays either death or survival benefit", {
  simulations <- tibble::tibble(
    Kx = c(0, 1, 3, 5)
  )

  result <- mc_insurance(
    simulations,
    i = 0,
    benefit = 5000,
    type = "endowment",
    n = 3,
    timing = "end_of_year"
  )

  expect_equal(result$benefit_indicator, rep(1, 4))
  expect_equal(result$benefit_time, c(1, 2, 3, 3))
  expect_equal(result$pv_benefit, rep(5000, 4))
})


test_that("mc_insurance moment-of-death benefits use complete lifetimes", {
  simulations <- tibble::tibble(
    Kx = c(0, 1, 3),
    Tx = c(0.4, 1.8, 3.2)
  )

  result <- mc_insurance(
    simulations,
    i = 0.05,
    benefit = 100,
    type = "term",
    n = 2.5,
    timing = "moment_of_death"
  )

  v <- 1 / 1.05

  expect_equal(result$benefit_indicator, c(1, 1, 0))
  expect_equal(
    result$pv_benefit,
    c(100 * v^0.4, 100 * v^1.8, 0),
    tolerance = 1e-12
  )
})


test_that("mc_insurance rejects fractional discrete terms and deferrals", {
  simulations <- tibble::tibble(
    Kx = c(1, 2),
    Tx = c(1.4, 2.6)
  )

  expect_error(
    mc_insurance(
      simulations,
      i = 0.05,
      type = "term",
      n = 2.5,
      timing = "end_of_year"
    ),
    "`n` must be a whole number"
  )

  expect_error(
    mc_insurance(
      simulations,
      i = 0.05,
      type = "deferred",
      h = 1.5,
      timing = "end_of_year"
    ),
    "`h` must be a whole number"
  )

  expect_no_error(
    mc_insurance(
      simulations,
      i = 0.05,
      type = "deferred_term",
      h = 1.5,
      n = 1.25,
      timing = "moment_of_death"
    )
  )
})


test_that("mc_insurance rejects noninteger interest conversion frequency", {
  simulations <- tibble::tibble(
    Kx = 2
  )

  expect_error(
    mc_insurance(
      simulations,
      i = 0.06,
      i_type = "nominal_interest",
      m = 1.5,
      type = "whole"
    ),
    "`m`"
  )
})


test_that("mc_insurance supports multiple-life status columns", {
  simulations <- tibble::tibble(
    K_status = c(1, 4),
    T_status = c(1.2, 4.7)
  )

  result <- mc_insurance(
    simulations,
    i = 0,
    benefit = 100000,
    type = "term",
    n = 3,
    timing = "moment_of_death",
    col_K = "K_status",
    col_T = "T_status"
  )

  expect_equal(result$benefit_indicator, c(1, 0))
  expect_equal(result$pv_benefit, c(100000, 0))
})


test_that("mc_insurance preserves missing-lifetime information", {
  simulations <- tibble::tibble(
    Kx = c(1, NA_real_)
  )

  result <- mc_insurance(
    simulations,
    i = 0.05,
    benefit = 100,
    type = "whole",
    timing = "end_of_year"
  )

  expect_equal(result$benefit_indicator[[1]], 1)
  expect_true(is.na(result$benefit_indicator[[2]]))
  expect_true(is.na(result$benefit_time[[2]]))
  expect_true(is.na(result$pv_benefit[[2]]))
})


test_that("mc_insurance deprecated aliases use exact conflict detection", {
  simulations <- tibble::tibble(
    Kx = 2
  )

  expect_error(
    mc_insurance(
      simulations,
      i = 0.05,
      h = 0,
      deferral_years = 0,
      type = "whole"
    ),
    "Provide only one of `h`"
  )

  result <- mc_insurance(
    data = simulations,
    rate = 0.05,
    insurance = "whole_life",
    payment_timing = "end_of_year"
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(result$type, "whole")
})
