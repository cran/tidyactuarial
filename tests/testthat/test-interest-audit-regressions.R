test_that("small positive and negative rates retain their first-order equivalents", {
  for (rate in c(1e-20, -1e-20)) {
    for (type in c("nominal_interest", "nominal_discount", "force")) {
      # Each equivalent is rate + O(rate^2). The correction is below
      # double precision relative to rate; division prevents a loose
      # absolute tolerance from accepting an erroneous zero.
      out <- standardize_interest(i_type = type, i = rate, m = 12)
      expect_equal(out / rate, 1, tolerance = 1e-14)
    }
    equivalents <- interest_equivalents(i = rate, m = 12)
    nominal_discount <- equivalents$value[equivalents$family == "nominal_discount"]
    expect_equal(nominal_discount / rate, 1, tolerance = 1e-14)
  }
})

test_that("missing rate conventions propagate without selecting a conversion", {
  out <- standardize_interest(i_type = c("effective", NA_character_, "force"),
    i = c(.05, .03, 0), m = c(1, NA_real_, 1))
  expect_equal(out, c(.05, NA_real_, 0))
  expect_true(is.na(standardize_interest(i_type = NA_character_, i = -.5)))
})

test_that("long-horizon accumulation preserves a tiny rate's material effect", {
  for (rate in c(1e-20, -1e-20)) {
    horizon <- 1e20
    # (1 +/- 1/n)^n -> exp(+/-1); finite-n relative error is O(1/n).
    expected <- exp(sign(rate))
    expect_equal(accumulation_factor(t = horizon, i = rate), expected,
                 tolerance = 1e-14)
    expect_equal(discount_factor(t = horizon, i = rate), 1 / expected,
                 tolerance = 1e-14)
    expect_equal(discount_factor_spot(t = horizon, i = rate), 1 / expected,
                 tolerance = 1e-14)
  }
})

test_that("accumulation satisfies interval composition for a nonconstant force", {
  # delta(t)=.03+.004t has primitive .03t+.002t^2.
  delta <- function(t) .03 + .004 * t
  primitive <- function(t) .03 * t + .002 * t^2
  expected <- exp(primitive(7) - primitive(2))
  expect_equal(accumulation_factor(s = 2, t = 7, delta = delta), expected)
  a <- accumulation_factor(s = 2, t = 4, delta = delta)
  b <- accumulation_factor(s = 4, t = 7, delta = delta)
  expect_equal(a * b, expected)
  expect_equal(discount_factor(s = 2, t = 7, delta = delta) * expected, 1)
})
