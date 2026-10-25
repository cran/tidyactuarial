test_that("md_table reproduces the annual multiple-decrement identities", {
  qx_df <- tibble::tibble(
    x = 30:33,
    q_death = c(0.10, 0.20, 0.30, 1.00),
    q_disability = c(0.20, 0.10, 0.10, 0.00)
  )

  md <- md_table(
    qx_df,
    radix = 1000,
    close = TRUE
  )

  expect_equal(md$q_total, c(0.30, 0.30, 0.40, 1.00), tolerance = 1e-12)
  expect_equal(md$p_total, c(0.70, 0.70, 0.60, 0.00), tolerance = 1e-12)

  # l_{x+1} = l_x p_x^(tau)
  expect_equal(md$lx, c(1000, 700, 490, 294), tolerance = 1e-12)

  # d_x^(tau) = l_x q_x^(tau)
  expect_equal(md$d_total, md$lx * md$q_total, tolerance = 1e-12)

  # Cause-specific decrements add to the total decrements.
  expect_equal(
    md$d_death + md$d_disability,
    md$d_total,
    tolerance = 1e-12
  )

  expect_equal(
    md$d_death,
    c(100, 140, 147, 294),
    tolerance = 1e-12
  )

  expect_equal(
    md$d_disability,
    c(200, 70, 49, 0),
    tolerance = 1e-12
  )
})

test_that("md_table sorts ages and supports an explicit age/cause specification", {
  qx_df <- tibble::tibble(
    age = c(32, 30, 33, 31),
    q_death = c(0.30, 0.10, 1.00, 0.20),
    q_disability = c(0.10, 0.20, 0.00, 0.10),
    label = c("c", "a", "d", "b")
  )

  md <- md_table(
    qx_df,
    age_col = "age",
    cause_cols = c("q_death", "q_disability"),
    radix = 1000,
    close = TRUE
  )

  expect_equal(md$x, 30:33)
  expect_equal(md$lx, c(1000, 700, 490, 294), tolerance = 1e-12)
  expect_false("label" %in% names(md))
})

test_that("md_table permits an open terminal age when close is FALSE", {
  qx_df <- tibble::tibble(
    x = 40:42,
    q_death = c(0.10, 0.15, 0.20),
    q_withdrawal = c(0.05, 0.05, 0.10)
  )

  md <- md_table(
    qx_df,
    radix = 1000,
    close = FALSE
  )

  expect_equal(md$q_total, c(0.15, 0.20, 0.30), tolerance = 1e-12)
  expect_equal(md$lx, c(1000, 850, 680), tolerance = 1e-12)
})

test_that("md_table validates the basic actuarial structure", {
  expect_error(
    md_table(
      tibble::tibble(
        x = 30:32,
        q_a = c(0.60, 0.20, 1.00),
        q_b = c(0.50, 0.10, 0.00)
      ),
      close = TRUE
    ),
    "q_total"
  )

  expect_error(
    md_table(
      tibble::tibble(
        x = c(30, 32, 33),
        q_a = c(0.10, 0.10, 1.00)
      ),
      close = TRUE
    ),
    "consecutive"
  )

  expect_error(
    md_table(
      tibble::tibble(
        x = 30:32,
        q_a = c(0.10, 0.20, 0.30)
      ),
      close = TRUE
    ),
    "last age"
  )
})
