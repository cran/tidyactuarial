make_lt_insurance_xy <- function() {
  lifetable(
    x = 60:68,
    lx = c(
      100000, 96000, 91000, 85000, 77000,
      67000, 54000, 36000, 0
    ),
    close = TRUE,
    frac = "UDD"
  )
}

test_that("insurance_xy term follows the standard annuity-due identity", {
  lt <- make_lt_insurance_xy()
  i <- 0.05
  v <- 1 / (1 + i)
  d <- i / (1 + i)
  n <- 4

  for (status in c("joint", "last")) {
    ann_due <- annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = i,
      status = status,
      n = n,
      k = 1,
      timing = "due",
      frac = "UDD"
    )

    p_n <- t_pxy(
      lt,
      x = 60,
      y = 62,
      t = n,
      frac = "UDD",
      status = status
    )

    expected <- 1 - d * ann_due - v^n * p_n

    expect_equal(
      insurance_xy(
        lt = lt,
        x = 60,
        y = 62,
        i = i,
        type = "term",
        status = status,
        n = n,
        frac = "UDD"
      ),
      expected,
      tolerance = 1e-12
    )
  }
})

test_that("insurance_xy endowment follows one minus d times the temporary annuity-due", {
  lt <- make_lt_insurance_xy()
  i <- 0.05
  d <- i / (1 + i)
  n <- 4

  for (status in c("joint", "last")) {
    ann_due <- annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = i,
      status = status,
      n = n,
      k = 1,
      timing = "due",
      frac = "UDD"
    )

    expected <- 1 - d * ann_due

    expect_equal(
      insurance_xy(
        lt = lt,
        x = 60,
        y = 62,
        i = i,
        type = "endowment",
        status = status,
        n = n,
        frac = "UDD"
      ),
      expected,
      tolerance = 1e-12
    )
  }
})

test_that("insurance_xy scales linearly with benefit", {
  lt <- make_lt_insurance_xy()

  unit <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    n = 4,
    frac = "UDD"
  )

  scaled <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    n = 4,
    benefit = 250000,
    frac = "UDD"
  )

  expect_equal(scaled, 250000 * unit, tolerance = 1e-8)
})

test_that("first-death insurance is at least as valuable as second-death insurance", {
  lt <- make_lt_insurance_xy()

  first_death <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    n = 4,
    frac = "UDD"
  )

  second_death <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "last",
    n = 4,
    frac = "UDD"
  )

  expect_gte(first_death, second_death)
})

test_that("insurance_xy deferment equals survival-discount times insurance at advanced ages", {
  lt <- make_lt_insurance_xy()
  i <- 0.05
  v <- 1 / (1 + i)
  h <- 1
  n <- 3

  deferred <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = i,
    type = "term",
    status = "joint",
    n = n,
    h = h,
    frac = "UDD"
  )

  continuation <- insurance_xy(
    lt = lt,
    x = 61,
    y = 63,
    i = i,
    type = "term",
    status = "joint",
    n = n,
    h = 0,
    frac = "UDD"
  )

  defer_factor <- v^h * t_pxy(
    lt,
    x = 60,
    y = 62,
    t = h,
    frac = "UDD",
    status = "joint"
  )

  expect_equal(
    deferred,
    defer_factor * continuation,
    tolerance = 1e-12
  )
})

test_that("insurance_xy tidy output contains the scalar APV", {
  lt <- make_lt_insurance_xy()

  scalar <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    n = 4,
    frac = "UDD"
  )

  out <- insurance_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    type = "term",
    status = "joint",
    n = 4,
    frac = "UDD",
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)

  numeric_candidates <- names(out)[vapply(out, is.numeric, logical(1))]
  expect_true(any(vapply(
    out[numeric_candidates],
    function(z) any(abs(z - scalar) < 1e-12, na.rm = TRUE),
    logical(1)
  )))
})
