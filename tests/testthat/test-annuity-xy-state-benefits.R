make_lt_annuity_xy <- function() {
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

test_that("annuity_xy joint-life annual immediate equals the direct APV sum", {
  lt <- make_lt_annuity_xy()
  i <- 0.05
  v <- 1 / (1 + i)

  expected <- sum(vapply(
    1:4,
    function(t) {
      px <- t_px(lt, x = 60, t = t, frac = "UDD")
      py <- t_px(lt, x = 62, t = t, frac = "UDD")
      v^t * px * py
    },
    numeric(1)
  ))

  expect_equal(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = i,
      status = "joint",
      n = 4,
      k = 1,
      timing = "immediate",
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_xy last-survivor annual immediate equals the direct APV sum", {
  lt <- make_lt_annuity_xy()
  i <- 0.05
  v <- 1 / (1 + i)

  expected <- sum(vapply(
    1:4,
    function(t) {
      px <- t_px(lt, x = 60, t = t, frac = "UDD")
      py <- t_px(lt, x = 62, t = t, frac = "UDD")
      v^t * (px + py - px * py)
    },
    numeric(1)
  ))

  expect_equal(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = i,
      status = "last",
      n = 4,
      k = 1,
      timing = "immediate",
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_xy state-based benefits match their expected-state payment formula", {
  lt <- make_lt_annuity_xy()
  i <- 0.04
  v <- 1 / (1 + i)

  b <- list(
    both = 1,
    x_only = 0.50,
    y_only = 0.25
  )

  expected <- sum(vapply(
    1:4,
    function(t) {
      px <- t_px(lt, x = 60, t = t, frac = "UDD")
      py <- t_px(lt, x = 62, t = t, frac = "UDD")

      expected_payment <- (
        b$both * px * py +
        b$x_only * px * (1 - py) +
        b$y_only * py * (1 - px)
      )

      v^t * expected_payment
    },
    numeric(1)
  ))

  expect_equal(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = i,
      benefit = b,
      n = 4,
      k = 1,
      timing = "immediate",
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_xy joint value does not exceed last-survivor value", {
  lt <- make_lt_annuity_xy()

  joint <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 4,
    frac = "UDD"
  )

  last <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "last",
    n = 4,
    frac = "UDD"
  )

  expect_lte(joint, last)
})

test_that("annuity_xy tidy output preserves APV and assumptions", {
  lt <- make_lt_annuity_xy()

  scalar <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 4,
    frac = "UDD"
  )

  out <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 4,
    frac = "UDD",
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)

  numeric_candidates <- names(out)[vapply(out, is.numeric, logical(1))]
  expect_true(length(numeric_candidates) >= 1L)
  expect_true(any(vapply(
    out[numeric_candidates],
    function(z) any(abs(z - scalar) < 1e-12, na.rm = TRUE),
    logical(1)
  )))
})
