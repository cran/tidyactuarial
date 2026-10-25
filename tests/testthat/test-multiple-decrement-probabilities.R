make_md_for_probability_tests <- function() {
  md_table(
    tibble::tibble(
      x = 30:33,
      q_death = c(0.10, 0.20, 0.30, 1.00),
      q_disability = c(0.20, 0.10, 0.10, 0.00)
    ),
    radix = 1000,
    close = TRUE
  )
}

test_that("t_qxj matches the direct integer-horizon cause-specific formula", {
  md <- make_md_for_probability_tests()

  # q_30^(death) + p_30^(tau) q_31^(death)
  expected_death <- 0.10 + 0.70 * 0.20

  expect_equal(
    t_qxj(md, x = 30, t = 2, cause = "q_death"),
    expected_death,
    tolerance = 1e-12
  )

  expected_disability <- 0.20 + 0.70 * 0.10

  expect_equal(
    t_qxj(md, x = 30, t = 2, cause = "q_disability"),
    expected_disability,
    tolerance = 1e-12
  )

  # Sum over mutually exclusive causes equals probability of any decrement.
  total_two_year <- 1 - 0.70 * 0.70

  expect_equal(
    expected_death + expected_disability,
    total_two_year,
    tolerance = 1e-12
  )
})

test_that("t_qxj returns zero at duration zero", {
  md <- make_md_for_probability_tests()

  expect_equal(
    t_qxj(md, x = 30, t = 0, cause = "q_death"),
    0,
    tolerance = 1e-12
  )
})

test_that("t_qxj implements UDD, constant-force and Balducci fractional tails", {
  md <- make_md_for_probability_tests()
  s <- 0.5

  # At age 30: q_total = 0.30 and q_death = 0.10.
  expected_udd <- s * 0.10

  expect_equal(
    t_qxj(md, x = 30, t = s, cause = "q_death", frac = "UDD"),
    expected_udd,
    tolerance = 1e-12
  )

  w_death <- 0.10 / 0.30
  expected_cf <- w_death * (1 - 0.70^s)

  expect_equal(
    t_qxj(md, x = 30, t = s, cause = "q_death", frac = "CF"),
    expected_cf,
    tolerance = 1e-12
  )

  q_tau_s_balducci <- (s * 0.30) / (1 - (1 - s) * 0.30)
  expected_balducci <- w_death * q_tau_s_balducci

  expect_equal(
    t_qxj(md, x = 30, t = s, cause = "q_death", frac = "Balducci"),
    expected_balducci,
    tolerance = 1e-12
  )
})

test_that("t_qxj combines an integer part with a fractional UDD tail", {
  md <- make_md_for_probability_tests()

  # First-year death plus survival to age 31 times half-year UDD death.
  expected <- 0.10 + 0.70 * (0.5 * 0.20)

  expect_equal(
    t_qxj(
      md,
      x = 30,
      t = 1.5,
      cause = "q_death",
      frac = "UDD"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("t_qxj vectorizes and provides a tidy audit output", {
  md <- make_md_for_probability_tests()

  out <- t_qxj(
    md,
    x = c(30, 31),
    t = c(1, 2),
    cause = "q_death",
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("x", "t", "cause", "frac", "tqxj"))

  expected <- c(
    0.10,
    0.20 + 0.70 * 0.30
  )

  expect_equal(out$tqxj, expected, tolerance = 1e-12)
  expect_true(all(is.na(out$frac)))
})

test_that("t_qxj rejects an unqualified fractional horizon and unknown cause", {
  md <- make_md_for_probability_tests()

  expect_error(
    t_qxj(md, x = 30, t = 1.5, cause = "q_death"),
    "integer-valued"
  )

  expect_error(
    t_qxj(md, x = 30, t = 1, cause = "q_unknown"),
    "not found"
  )
})

test_that("lt_tau reproduces the total-decrement life table", {
  md <- make_md_for_probability_tests()

  lt <- lt_tau(
    md,
    radix = 1000,
    close = TRUE,
    frac = "UDD"
  )

  expect_s3_class(lt, "lifetable")
  expect_equal(lt$x, md$x)
  expect_equal(lt$lx, md$lx, tolerance = 1e-12)
  expect_equal(lt$qx, md$q_total, tolerance = 1e-12)
  expect_equal(lt$px, md$p_total, tolerance = 1e-12)

  # Two-year survival against all decrements.
  expect_equal(
    t_px(lt, x = 30, t = 2),
    0.70 * 0.70,
    tolerance = 1e-12
  )
})
