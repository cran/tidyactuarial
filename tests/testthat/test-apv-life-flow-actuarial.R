make_lt_apv_flow <- function() {
  tibble::tibble(
    x = 60:64,
    lx = c(1000, 800, 600, 300, 0)
  )
}

test_that("apv_life_flow matches direct single-life survival-discount valuation", {
  lt <- make_lt_apv_flow()
  t <- 0:3
  cf <- c(100, 100, 100, 100)
  i <- 0.05

  expected_survival <- c(1.0, 0.8, 0.6, 0.3)
  expected_discount <- (1 + i)^(-t)
  expected_pv <- cf * expected_survival * expected_discount

  out <- apv_life_flow(
    lt = lt,
    ages = 60,
    t = t,
    cf = cf,
    i = i
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$surv_prob, expected_survival, tolerance = 1e-12)
  expect_equal(out$discount, expected_discount, tolerance = 1e-12)
  expect_equal(out$expected_cf, cf * expected_survival, tolerance = 1e-12)
  expect_equal(out$pv, expected_pv, tolerance = 1e-12)
  expect_equal(out$pv_cum, cumsum(expected_pv), tolerance = 1e-12)
  expect_equal(attr(out, "apv"), sum(expected_pv), tolerance = 1e-12)
})

test_that("apv_life_flow uses UDD for fractional single-life survival", {
  lt <- make_lt_apv_flow()
  i <- 0.04

  # q60 = .2, so under UDD S(.5)=1-.5(.2)=.9.
  out <- apv_life_flow(
    lt = lt,
    ages = 60,
    t = 0.5,
    cf = 1000,
    i = i
  )

  expect_equal(out$surv_prob, 0.90, tolerance = 1e-12)
  expect_equal(
    out$pv,
    1000 * 0.90 * (1 + i)^(-0.5),
    tolerance = 1e-10
  )
})

test_that("apv_life_flow first-death status equals the product of individual survivals", {
  lt <- make_lt_apv_flow()
  t <- c(1, 2)

  # Life 60: S(1)=.8, S(2)=.6
  # Life 61: S(1)=600/800=.75, S(2)=300/800=.375
  expected <- c(
    0.8 * 0.75,
    0.6 * 0.375
  )

  out <- apv_life_flow(
    lt = lt,
    ages = c(60, 61),
    t = t,
    cf = c(1, 1),
    i = 0,
    status = "first"
  )

  expect_equal(out$surv_prob, expected, tolerance = 1e-12)
  expect_equal(attr(out, "apv"), sum(expected), tolerance = 1e-12)
})

test_that("apv_life_flow last-survivor status uses inclusion-exclusion", {
  lt <- make_lt_apv_flow()
  t <- c(1, 2)

  p60 <- c(0.8, 0.6)
  p61 <- c(0.75, 0.375)

  expected <- 1 - (1 - p60) * (1 - p61)

  out <- apv_life_flow(
    lt = lt,
    ages = c(60, 61),
    t = t,
    cf = c(1, 1),
    i = 0,
    status = "last"
  )

  expect_equal(out$surv_prob, expected, tolerance = 1e-12)
})

test_that("apv_life_flow reversionary status blends joint and last-survivor states", {
  lt <- make_lt_apv_flow()
  t <- c(1, 2)
  alpha <- 0.40

  p60 <- c(0.8, 0.6)
  p61 <- c(0.75, 0.375)

  p_all <- p60 * p61
  p_any <- 1 - (1 - p60) * (1 - p61)
  expected <- p_all + alpha * (p_any - p_all)

  out <- apv_life_flow(
    lt = lt,
    ages = c(60, 61),
    t = t,
    cf = c(1, 1),
    i = 0,
    status = "reversionary",
    alpha = alpha
  )

  expect_equal(out$surv_prob, expected, tolerance = 1e-12)
})

test_that("apv_life_flow respects equivalent nominal interest inputs", {
  lt <- make_lt_apv_flow()
  t <- 1:3
  cf <- c(100, 200, 300)

  j12 <- 0.06
  i_eff <- (1 + j12 / 12)^12 - 1

  nominal <- apv_life_flow(
    lt = lt,
    ages = 60,
    t = t,
    cf = cf,
    i = j12,
    i_type = "nominal_interest",
    m = 12
  )

  effective <- apv_life_flow(
    lt = lt,
    ages = 60,
    t = t,
    cf = cf,
    i = i_eff,
    i_type = "effective"
  )

  expect_equal(nominal$discount, effective$discount, tolerance = 1e-12)
  expect_equal(nominal$pv, effective$pv, tolerance = 1e-10)
  expect_equal(attr(nominal, "apv"), attr(effective, "apv"), tolerance = 1e-10)
})

test_that("apv_life_flow supports payment dates and an explicit valuation date", {
  lt <- make_lt_apv_flow()

  date0 <- as.Date("2026-01-01")
  dates <- as.Date(c("2026-01-01", "2027-01-01"))
  times <- as.numeric(dates - date0) / 365.25

  out <- apv_life_flow(
    lt = lt,
    ages = 60,
    date = dates,
    date0 = date0,
    cf = c(100, 100),
    i = 0.05
  )

  expect_equal(out$date, dates)
  expect_equal(out$t, times, tolerance = 1e-12)

  expected_surv_second <- 1 - times[2] * 0.20
  expect_equal(
    out$surv_prob,
    c(1, expected_surv_second),
    tolerance = 1e-12
  )
})

test_that("apv_life_flow attaches a cumulative APV plot on request", {
  lt <- make_lt_apv_flow()

  out <- apv_life_flow(
    lt = lt,
    ages = 60,
    t = 0:3,
    cf = rep(100, 4),
    i = 0.05,
    plot = TRUE
  )

  p <- attr(out, "plot")

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 2L)
})

test_that("apv_life_flow validates status-specific and timing inputs", {
  lt <- make_lt_apv_flow()

  expect_error(
    apv_life_flow(
      lt = lt,
      ages = c(60, 61),
      t = 1,
      cf = 100,
      i = 0.05,
      status = "reversionary"
    ),
    "`alpha`"
  )

  expect_error(
    apv_life_flow(
      lt = lt,
      ages = 60,
      t = 1,
      date = as.Date("2027-01-01"),
      cf = 100,
      i = 0.05
    ),
    "only one"
  )

  expect_error(
    apv_life_flow(
      lt = lt,
      ages = 60,
      t = c(1, 2),
      cf = 100,
      i = 0.05
    ),
    "same length"
  )
})
