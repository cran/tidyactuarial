test_that("reserve_xy supports two different life tables in one call", {
  lt_x <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.011 * (0:60)^1.34))
  )

  lt_y <- data.frame(
    x = 40:105,
    lx = round(100000 * exp(-0.014 * (0:65)^1.37))
  )

  lt_x$lx[nrow(lt_x)] <- 0
  lt_y$lx[nrow(lt_y)] <- 0

  result <- reserve_xy(
    lt = list(lt_x, lt_y),
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 12,
    n_prem = 10,
    frac = "UDD",
    t = c(0, 5, 10, 15, 20),
    output = "summary"
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 5L)
  expect_true(all(is.finite(result$reserve)))
  expect_equal(
    result$premium_annualized,
    12 * result$premium_per_payment,
    tolerance = 1e-12
  )
})


test_that("two-table direct and piped reserve_xy calculations agree", {
  lt_x <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.011 * (0:60)^1.34))
  )

  lt_y <- data.frame(
    x = 40:105,
    lx = round(100000 * exp(-0.014 * (0:65)^1.37))
  )

  lt_x$lx[nrow(lt_x)] <- 0
  lt_y$lx[nrow(lt_y)] <- 0

  direct <- reserve_xy(
    lt = list(lt_x, lt_y),
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "last",
    benefit = 100000,
    n = 20,
    k = 4,
    n_prem = 10,
    frac = "UDD",
    t = c(0, 5, 10, 15, 20),
    output = "value"
  )

  piped <- life_contract(
    lt = list(lt_x, lt_y),
    lives = "last_survivor",
    x = 60,
    y = 62,
    i = 0.05
  ) |>
    add_insurance(
      type = "term",
      benefit = 100000,
      n = 20,
      frac = "UDD"
    ) |>
    add_premium_schedule(
      k = 4,
      n_prem = 10,
      frac = "UDD"
    ) |>
    reserve_xy(
      t = c(0, 5, 10, 15, 20),
      output = "value"
    )

  expect_equal(piped, direct, tolerance = 1e-10)
})


test_that("reserve_xy applies each life table to its corresponding life", {
  lt_x <- data.frame(
    x = 40:100,
    lx = round(100000 * exp(-0.010 * (0:60)^1.30))
  )

  lt_y_low_mortality <- data.frame(
    x = 40:105,
    lx = round(100000 * exp(-0.008 * (0:65)^1.28))
  )

  lt_y_high_mortality <- data.frame(
    x = 40:105,
    lx = round(100000 * exp(-0.022 * (0:65)^1.45))
  )

  lt_x$lx[nrow(lt_x)] <- 0
  lt_y_low_mortality$lx[nrow(lt_y_low_mortality)] <- 0
  lt_y_high_mortality$lx[nrow(lt_y_high_mortality)] <- 0

  reserve_low <- reserve_xy(
    lt = list(lt_x, lt_y_low_mortality),
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = 5,
    output = "value"
  )

  reserve_high <- reserve_xy(
    lt = list(lt_x, lt_y_high_mortality),
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    benefit = 100000,
    n = 20,
    k = 1,
    n_prem = 10,
    t = 5,
    output = "value"
  )

  expect_false(isTRUE(all.equal(reserve_low, reserve_high)))
})
