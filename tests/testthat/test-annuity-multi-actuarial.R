make_lt_annuity_multi <- function() {
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

test_that("annuity_multi joint-life immediate equals a direct discounted survival sum", {
  lt <- make_lt_annuity_multi()
  i <- 0.05
  v <- 1 / (1 + i)
  n <- 4

  expected <- sum(vapply(
    seq_len(n),
    function(t) {
      v^t *
        t_px(lt, x = 60, t = t, frac = "UDD") *
        t_px(lt, x = 62, t = t, frac = "UDD")
    },
    numeric(1)
  ))

  expect_equal(
    annuity_multi(
      lt = lt,
      ages = c(60, 62),
      i = i,
      n = n,
      annuity = "cohort",
      cohort = "first",
      timing = "immediate"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_multi last-survivor immediate equals inclusion-exclusion at each payment time", {
  lt <- make_lt_annuity_multi()
  i <- 0.04
  v <- 1 / (1 + i)
  n <- 4

  expected <- sum(vapply(
    seq_len(n),
    function(t) {
      px <- t_px(lt, x = 60, t = t, frac = "UDD")
      py <- t_px(lt, x = 62, t = t, frac = "UDD")
      v^t * (px + py - px * py)
    },
    numeric(1)
  ))

  expect_equal(
    annuity_multi(
      lt = lt,
      ages = c(60, 62),
      i = i,
      n = n,
      annuity = "cohort",
      cohort = "last",
      timing = "immediate"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("annuity_multi supports three-life joint and last-survivor statuses", {
  lt <- make_lt_annuity_multi()
  i <- 0.03
  v <- 1 / (1 + i)
  n <- 3
  ages <- c(60, 61, 62)

  expected_joint <- sum(vapply(
    seq_len(n),
    function(t) {
      p <- vapply(
        ages,
        function(age) t_px(lt, x = age, t = t, frac = "UDD"),
        numeric(1)
      )
      v^t * prod(p)
    },
    numeric(1)
  ))

  expected_last <- sum(vapply(
    seq_len(n),
    function(t) {
      p <- vapply(
        ages,
        function(age) t_px(lt, x = age, t = t, frac = "UDD"),
        numeric(1)
      )
      v^t * (1 - prod(1 - p))
    },
    numeric(1)
  ))

  joint <- annuity_multi(
    lt = lt,
    ages = ages,
    i = i,
    n = n,
    annuity = "cohort",
    cohort = "first",
    timing = "immediate"
  )

  last <- annuity_multi(
    lt = lt,
    ages = ages,
    i = i,
    n = n,
    annuity = "cohort",
    cohort = "last",
    timing = "immediate"
  )

  expect_equal(joint, expected_joint, tolerance = 1e-12)
  expect_equal(last, expected_last, tolerance = 1e-12)
  expect_gte(last, joint)
})

test_that("annuity_multi reversionary value lies between joint and last survivor", {
  lt <- make_lt_annuity_multi()
  i <- 0.05
  n <- 4
  alpha <- 0.40

  joint <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = i,
    n = n,
    annuity = "cohort",
    cohort = "first",
    timing = "immediate"
  )

  last <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = i,
    n = n,
    annuity = "cohort",
    cohort = "last",
    timing = "immediate"
  )

  reversionary <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = i,
    n = n,
    annuity = "reversionary",
    alpha = alpha,
    timing = "immediate"
  )

  expect_equal(
    reversionary,
    joint + alpha * (last - joint),
    tolerance = 1e-12
  )

  expect_gte(reversionary, joint)
  expect_lte(reversionary, last)
})

test_that("annuity_multi alpha endpoints reproduce joint and last survivor", {
  lt <- make_lt_annuity_multi()

  joint <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = 0.05,
    n = 4,
    annuity = "cohort",
    cohort = "first",
    timing = "immediate"
  )

  last <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = 0.05,
    n = 4,
    annuity = "cohort",
    cohort = "last",
    timing = "immediate"
  )

  expect_equal(
    annuity_multi(
      lt = lt,
      ages = c(60, 62),
      i = 0.05,
      n = 4,
      annuity = "reversionary",
      alpha = 0,
      timing = "immediate"
    ),
    joint,
    tolerance = 1e-12
  )

  expect_equal(
    annuity_multi(
      lt = lt,
      ages = c(60, 62),
      i = 0.05,
      n = 4,
      annuity = "reversionary",
      alpha = 1,
      timing = "immediate"
    ),
    last,
    tolerance = 1e-12
  )
})

test_that("annuity_multi due timing adds the time-zero status payment", {
  lt <- make_lt_annuity_multi()
  i <- 0.05
  n <- 4

  immediate <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = i,
    n = n,
    annuity = "cohort",
    cohort = "first",
    timing = "immediate"
  )

  due <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = i,
    n = n,
    annuity = "cohort",
    cohort = "first",
    timing = "due"
  )

  # For a finite n-year annuity, the immediate and due payment sets are shifted:
  # due pays at t=0,...,n-1; immediate at t=1,...,n.
  expected_due <- sum(vapply(
    0:(n - 1),
    function(t) {
      (1 + i)^(-t) *
        t_px(lt, x = 60, t = t, frac = "UDD") *
        t_px(lt, x = 62, t = t, frac = "UDD")
    },
    numeric(1)
  ))

  expect_equal(due, expected_due, tolerance = 1e-12)
  expect_gt(due, immediate)
})

test_that("annuity_multi matches annuity_xy for two-life annual cohort annuities", {
  lt <- make_lt_annuity_multi()

  multi_joint <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = 0.05,
    n = 4,
    annuity = "cohort",
    cohort = "first",
    timing = "immediate"
  )

  xy_joint <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "joint",
    n = 4,
    k = 1,
    timing = "immediate",
    frac = "UDD"
  )

  multi_last <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = 0.05,
    n = 4,
    annuity = "cohort",
    cohort = "last",
    timing = "immediate"
  )

  xy_last <- annuity_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    status = "last",
    n = 4,
    k = 1,
    timing = "immediate",
    frac = "UDD"
  )

  expect_equal(multi_joint, xy_joint, tolerance = 1e-12)
  expect_equal(multi_last, xy_last, tolerance = 1e-12)
})

test_that("annuity_multi respects equivalent interest-rate specifications", {
  lt <- make_lt_annuity_multi()

  j12 <- 0.06
  i_eff <- (1 + j12 / 12)^12 - 1

  nominal <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = j12,
    i_type = "nominal_interest",
    m = 12,
    n = 4,
    annuity = "cohort",
    cohort = "first"
  )

  effective <- annuity_multi(
    lt = lt,
    ages = c(60, 62),
    i = i_eff,
    i_type = "effective",
    n = 4,
    annuity = "cohort",
    cohort = "first"
  )

  expect_equal(nominal, effective, tolerance = 1e-12)
})
