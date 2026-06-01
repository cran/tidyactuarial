test_that("annuity_x computes annual annuity-due under compact actuarial notation", {
  lt <- data.frame(
    x = 60:65,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000)
  )

  ax <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.06,
    n = 3,
    timing = "due"
  )

  v <- 1 / 1.06
  expected <- sum(
    v^(0:2) * lt$lx[match(60 + 0:2, lt$x)] / lt$lx[lt$x == 60]
  )

  expect_type(ax, "double")
  expect_length(ax, 1)
  expect_equal(ax, expected, tolerance = 1e-10)
})

test_that("annuity_x computes annual annuity-immediate under compact actuarial notation", {
  lt <- data.frame(
    x = 60:65,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000)
  )

  ax <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.06,
    n = 3,
    timing = "immediate"
  )

  v <- 1 / 1.06
  expected <- sum(
    v^(1:3) * lt$lx[match(60 + 1:3, lt$x)] / lt$lx[lt$x == 60]
  )

  expect_type(ax, "double")
  expect_length(ax, 1)
  expect_equal(ax, expected, tolerance = 1e-10)
})

test_that("annuity_x handles deferment with h", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  deferred <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 3,
    h = 2,
    timing = "due"
  )

  v <- 1 / 1.05

  deferment_factor <- v^2 * lt$lx[lt$x == 62] / lt$lx[lt$x == 60]

  annuity_at_62 <- sum(
    v^(0:2) * lt$lx[match(62 + 0:2, lt$x)] / lt$lx[lt$x == 62]
  )

  expected <- deferment_factor * annuity_at_62

  expect_equal(deferred, expected, tolerance = 1e-10)
})

test_that("annuity_x supports k-thly payments with fractional assumption", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  ax <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 3,
    k = 12,
    timing = "due",
    frac = "UDD",
    woolhouse = "none"
  )

  expect_type(ax, "double")
  expect_length(ax, 1)
  expect_true(is.finite(ax))
  expect_gt(ax, 0)
})

test_that("annuity_x supports nominal rates through i_type and m", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  ax_nominal <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.06,
    i_type = "nominal_interest",
    m = 12,
    n = 3,
    timing = "due"
  )

  i_eff <- (1 + 0.06 / 12)^12 - 1

  ax_effective <- annuity_x(
    lt = lt,
    x = 60,
    i = i_eff,
    i_type = "effective",
    m = 1,
    n = 3,
    timing = "due"
  )

  expect_equal(ax_nominal, ax_effective, tolerance = 1e-10)
})

test_that("annuity_x returns tidy output when tidy is TRUE", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  out <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 3,
    h = 1,
    k = 1,
    timing = "due",
    tidy = TRUE
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c(
    "x",
    "i",
    "i_type",
    "m",
    "n",
    "h",
    "k",
    "timing",
    "deferment_factor",
    "pure_endowment_factor",
    "apv"
  ) %in% names(out)))
  expect_equal(nrow(out), 1)
  expect_true(is.finite(out$apv))
})

test_that("annuity_x works with a single-life life_contract", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  contract <- life_contract(
    lt = lt,
    lives = "single",
    x = 60,
    i = 0.05
  )

  via_contract <- contract |>
    annuity_x(
      n = 3,
      timing = "due"
    )

  direct <- annuity_x(
    lt = lt,
    x = 60,
    i = 0.05,
    n = 3,
    timing = "due"
  )

  expect_equal(via_contract, direct, tolerance = 1e-10)
})

test_that("annuity_x rejects old argument names", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  expect_error(
    annuity_x(
      mortality_table = lt,
      age = 60,
      rate = 0.05,
      term_years = 3
    )
  )
})

test_that("annuity_x validates compact actuarial arguments", {
  lt <- data.frame(
    x = 60:70,
    lx = seq(100000, 80000, length.out = 11)
  )

  expect_error(
    annuity_x(lt = lt, x = 60, i = 0.05, h = -1),
    "`h` must be a single nonnegative integer"
  )

  expect_error(
    annuity_x(lt = lt, x = 60, i = 0.05, k = 0),
    "`k` must be a single positive integer"
  )

  expect_error(
    annuity_x(lt = lt, x = 60, i = 0.05, i_type = "bad"),
    "`i_type` must be one of"
  )

  expect_error(
    annuity_x(lt = lt, x = 60, i = 0.05, n = 20),
    "`n` exceeds the horizon"
  )
})
