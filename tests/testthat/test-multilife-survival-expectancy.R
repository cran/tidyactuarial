make_lt_multilife <- function() {
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

test_that("t_pxy joint and last-survivor probabilities satisfy independence identities", {
  lt <- make_lt_multilife()

  px <- t_px(lt, x = 60, t = 2.5, frac = "UDD")
  py <- t_px(lt, x = 62, t = 2.5, frac = "UDD")

  joint <- t_pxy(
    lt,
    x = 60,
    y = 62,
    t = 2.5,
    frac = "UDD",
    status = "joint"
  )

  last <- t_pxy(
    lt,
    x = 60,
    y = 62,
    t = 2.5,
    frac = "UDD",
    status = "last"
  )

  expect_equal(joint, px * py, tolerance = 1e-12)
  expect_equal(last, px + py - px * py, tolerance = 1e-12)
})

test_that("t_pxy accepts two different mortality tables", {
  lt_x <- lifetable(
    x = 60:68,
    lx = c(100000, 96000, 91000, 85000, 77000, 67000, 54000, 36000, 0),
    close = TRUE,
    frac = "UDD"
  )

  lt_y <- lifetable(
    x = 60:68,
    lx = c(100000, 97500, 94000, 89500, 83500, 75500, 65000, 50000, 0),
    close = TRUE,
    frac = "UDD"
  )

  px <- t_px(lt_x, x = 60, t = 3, frac = "UDD")
  py <- t_px(lt_y, x = 62, t = 3, frac = "UDD")

  expect_equal(
    t_pxy(
      list(lt_x, lt_y),
      x = 60,
      y = 62,
      t = 3,
      frac = "UDD",
      status = "joint"
    ),
    px * py,
    tolerance = 1e-12
  )
})

test_that("e_xy curtate joint-life expectation equals the direct survival sum", {
  lt <- make_lt_multilife()
  horizon <- 4

  expected <- sum(vapply(
    seq_len(horizon),
    function(k) {
      t_pxy(
        lt,
        x = 60,
        y = 62,
        t = k,
        frac = "UDD",
        status = "joint"
      )
    },
    numeric(1)
  ))

  expect_equal(
    e_xy(
      lt,
      x = 60,
      y = 62,
      t = horizon,
      type = "curtate",
      cohort = "first"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("e_xy curtate last-survivor expectation equals the direct survival sum", {
  lt <- make_lt_multilife()
  horizon <- 4

  expected <- sum(vapply(
    seq_len(horizon),
    function(k) {
      t_pxy(
        lt,
        x = 60,
        y = 62,
        t = k,
        frac = "UDD",
        status = "last"
      )
    },
    numeric(1)
  ))

  expect_equal(
    e_xy(
      lt,
      x = 60,
      y = 62,
      t = horizon,
      type = "curtate",
      cohort = "last"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("e_xy obeys the last-survivor inclusion-exclusion identity", {
  lt <- make_lt_multilife()
  horizon <- 4

  ex <- e_x(
    lt,
    x = 60,
    t = horizon,
    type = "curtate"
  )

  ey <- e_x(
    lt,
    x = 62,
    t = horizon,
    type = "curtate"
  )

  e_joint <- e_xy(
    lt,
    x = 60,
    y = 62,
    t = horizon,
    type = "curtate",
    cohort = "first"
  )

  e_last <- e_xy(
    lt,
    x = 60,
    y = 62,
    t = horizon,
    type = "curtate",
    cohort = "last"
  )

  expect_equal(
    e_last,
    ex + ey - e_joint,
    tolerance = 1e-12
  )
})

test_that("e_xy returns a tidy row containing the scalar actuarial result", {
  lt <- make_lt_multilife()

  scalar <- e_xy(
    lt,
    x = 60,
    y = 62,
    t = 3,
    type = "curtate",
    cohort = "first"
  )

  out <- e_xy(
    lt,
    x = 60,
    y = 62,
    t = 3,
    type = "curtate",
    cohort = "first",
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)

  # The public documentation promises a tibble for tidy = TRUE, but it does
  # not prescribe the exact name of the result column. Therefore the test
  # verifies the contract without imposing an undocumented column name.
  numeric_cols <- names(out)[vapply(out, is.numeric, logical(1))]

  expect_true(length(numeric_cols) >= 1L)
  expect_true(any(vapply(
    out[numeric_cols],
    function(z) any(abs(z - scalar) < 1e-12, na.rm = TRUE),
    logical(1)
  )))
})
