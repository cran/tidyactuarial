test_that("annuity_xy exact k-thly valuation uses payments of 1 / k", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  i <- 0.05
  k <- 12L
  n <- 2

  times <- (0:(k * n - 1L)) / k

  expected <- sum(
    (1 / k) *
      (1 + i)^(-times) *
      vapply(
        times,
        function(tt) {
          px <- t_px(
            lt = lt,
            x = 60,
            t = tt,
            frac = "UDD",
            tidy = FALSE,
            check = FALSE
          )

          py <- t_px(
            lt = lt,
            x = 62,
            t = tt,
            frac = "UDD",
            tidy = FALSE,
            check = FALSE
          )

          px * py
        },
        numeric(1L)
      )
  )

  observed <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = i,
    status = "joint",
    n = n,
    k = k,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("annuity_xy accepts fractional exact terms aligned with k", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  result <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 2.5,
    k = 12,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  expect_true(is.finite(result))
  expect_gt(result, 0)
})


test_that("annuity_xy rejects exact terms not aligned with k", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  expect_error(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = 0.05,
      status = "joint",
      n = 2.3,
      k = 4,
      woolhouse = "none"
    ),
    "`n \\* k` must be an integer"
  )
})


test_that("annuity_xy keeps Woolhouse terms integer-valued", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  expect_error(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = 0.05,
      status = "joint",
      n = 2.5,
      k = 12,
      woolhouse = "first"
    ),
    "Woolhouse approximations require `n` to be an integer"
  )
})


test_that("annuity_xy preserves the deferred due-immediate identity", {
  lt <- data.frame(
    x = 60:75,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000, 51000,
      42000, 33000, 24000, 15000
    )
  )

  i <- 0.05
  h <- 2
  n <- 2.5
  k <- 12

  due <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = i,
    status = "joint",
    n = n,
    h = h,
    k = k,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  immediate <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = i,
    status = "joint",
    n = n,
    h = h,
    k = k,
    timing = "immediate",
    woolhouse = "none",
    frac = "UDD"
  )

  status_probability <- function(tt) {
    px <- t_px(
      lt = lt,
      x = 60,
      t = tt,
      frac = "UDD",
      tidy = FALSE,
      check = FALSE
    )

    py <- t_px(
      lt = lt,
      x = 62,
      t = tt,
      frac = "UDD",
      tidy = FALSE,
      check = FALSE
    )

    px * py
  }

  endpoint_start <- (1 + i)^(-h) * status_probability(h)
  endpoint_end <- (1 + i)^(-(h + n)) * status_probability(h + n)

  expect_equal(
    immediate,
    due - (1 / k) * (endpoint_start - endpoint_end),
    tolerance = 1e-12
  )
})


test_that("annuity_xy first-order Woolhouse scales deferred endpoints", {
  lt <- data.frame(
    x = 60:75,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000, 51000,
      42000, 33000, 24000, 15000
    )
  )

  i <- 0.05
  h <- 2
  n <- 4
  k <- 12

  annual_due <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = i,
    status = "joint",
    n = n,
    h = h,
    k = 1,
    timing = "due",
    woolhouse = "none",
    frac = "UDD"
  )

  status_probability <- function(tt) {
    px <- t_px(
      lt = lt,
      x = 60,
      t = tt,
      frac = "UDD",
      tidy = FALSE,
      check = FALSE
    )

    py <- t_px(
      lt = lt,
      x = 62,
      t = tt,
      frac = "UDD",
      tidy = FALSE,
      check = FALSE
    )

    px * py
  }

  endpoint_start <- (1 + i)^(-h) * status_probability(h)
  endpoint_end <- (1 + i)^(-(h + n)) * status_probability(h + n)

  expected <- annual_due -
    (k - 1) / (2 * k) * (endpoint_start - endpoint_end)

  observed <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = i,
    status = "joint",
    n = n,
    h = h,
    k = k,
    timing = "due",
    woolhouse = "first",
    frac = "UDD"
  )

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("annuity_xy restricts Woolhouse for state-based benefits", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  expect_error(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = 0.05,
      benefit = list(
        both = 0,
        x_only = 1,
        y_only = 0
      ),
      n = 4,
      k = 12,
      woolhouse = "first"
    ),
    "Woolhouse approximations are supported only"
  )
})


test_that("annuity_xy contract and direct calculations agree", {
  lt <- data.frame(
    x = 60:70,
    lx = c(
      100000, 99000, 97500, 95500, 93000, 90000,
      86000, 81000, 75000, 68000, 60000
    )
  )

  direct <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 4,
    k = 12,
    timing = "due",
    frac = "UDD"
  )

  piped <- life_contract(
    lt = lt,
    lives = "joint",
    x = 60,
    y = 62,
    i = 0.05
  ) |>
    annuity_xy(
      n = 4,
      k = 12,
      timing = "due",
      frac = "UDD"
    )

  expect_equal(piped, direct, tolerance = 1e-12)
})


test_that("annuity_xy rejects terms beyond the available status horizon", {
  lt <- data.frame(
    x = 60:66,
    lx = c(100000, 99000, 97500, 95500, 93000, 90000, 0)
  )

  expect_error(
    annuity_xy(
      lt = lt,
      x = 60,
      y = 62,
      i = 0.05,
      status = "joint",
      n = 10,
      h = 1,
      k = 1
    ),
    "exceeds the horizon supported"
  )
})
