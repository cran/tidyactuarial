test_that("complete last-survivor expectation preserves the initial conditioning", {
  lt <- data.frame(x = 60:64, lx = c(100, 80, 50, 20, 0))
  # Independent UDD reference: each marginal survival is piecewise linear.
  sx <- approxfun(0:4, lt$lx / 100, yright = 0, rule = 2)
  sy <- approxfun(0:3, lt$lx[2:5] / 80, yright = 0, rule = 2)
  for (horizon in c(0, 0.4, 1, 2.5, 4, 8)) {
    end <- min(horizon, 4)
    breaks <- sort(unique(c(0, seq_len(floor(end)), end)))
    breaks <- breaks[breaks <= end]
    joint <- last <- 0
    if (end > 0) {
      for (j in seq_len(length(breaks) - 1L)) {
        a <- breaks[j]; b <- breaks[j + 1L]
        joint <- joint + integrate(function(t) sx(t) * sy(t), a, b)$value
        last <- last + integrate(function(t) sx(t) + sy(t) - sx(t) * sy(t), a, b)$value
      }
    }
    expect_equal(e_xy(lt, 60, 61, t = horizon, type = "complete",
                       cohort = "first"), joint, tolerance = 1e-10)
    expect_equal(e_xy(lt, 60, 61, t = horizon, type = "complete",
                       cohort = "last"), last, tolerance = 1e-10)
  }
})

test_that("complete two-life expectations satisfy inclusion-exclusion for all fractions", {
  lt <- lifetable(60:65, lx = c(100, 90, 70, 40, 20, 10), close = FALSE)
  for (frac in c("UDD", "CF", "CML", "Balducci")) {
    for (horizon in c(2.5, 5)) {
      joint <- e_xy(lt, 60, 62, t = horizon, type = "complete", frac = frac)
      last <- e_xy(lt, 60, 62, t = horizon, type = "complete", frac = frac,
                    cohort = "last")
      reference <- e_x(lt, 60, t = horizon, type = "complete", frac = frac) +
        e_x(lt, 62, t = horizon, type = "complete", frac = frac)
      expect_equal(joint + last, reference, tolerance = 1e-9)
    }
  }
})

test_that("integer and extinct survival do not require a next-year probability", {
  lt <- data.frame(x = 60:64, lx = c(100, 80, 50, 0, 0))
  expect_equal(t_px(lt, 60, 4), 0)
  expect_equal(t_px(lt, 60, 3.5), 0)
  expect_equal(t_px(data.frame(x = 60:61, lx = c(100, 80)), 60, 1), 0.8)
  expect_equal(e_xy(lt, 60, 61, cohort = "first"), 0.5)
  expect_equal(e_xy(lt, 60, 61, cohort = "last"), 1.425)
})

test_that("e_xy handles invalid initial lives before returning a zero horizon", {
  lt <- data.frame(x = 60:64, lx = c(100, 80, 50, 20, 0))
  expect_true(is.na(e_xy(lt, 60, 64)))
  expect_true(is.na(e_xy(lt, 60, 65)))
  expect_true(is.na(e_xy(lt, 59, 61, t = 0)))
})

test_that("complete CF and Balducci results agree with analytic integrals", {
  mu <- 0.2
  lt_cf <- data.frame(x = 60:64, lx = 100 * exp(-mu * (0:4)))
  # Different horizons: joint ends after 3 years, last survivor after 4.
  joint_cf <- -expm1(-2 * mu * 3) / (2 * mu)
  last_cf <- -expm1(-mu * 4) / mu - expm1(-mu * 3) / mu - joint_cf
  expect_equal(e_xy(lt_cf, 60, 61, type = "complete", frac = "CF"),
               joint_cf, tolerance = 1e-10)
  expect_equal(e_xy(lt_cf, 60, 61, type = "complete", frac = "CML",
                    cohort = "last"), last_cf, tolerance = 1e-10)

  # Reciprocal lx is exactly Balducci: S_60(t) = 1/(1+t),
  # S_61(t) = 2/(2+t), until the respective table horizons.
  lt_bal <- data.frame(x = 60:64, lx = 100 / (1:5))
  joint_bal <- 2 * (log1p(3) - log1p(3 / 2))
  last_bal <- log1p(4) + 2 * log1p(3 / 2) - joint_bal
  expect_equal(e_xy(lt_bal, 60, 61, type = "complete", frac = "Balducci"),
               joint_bal, tolerance = 1e-10)
  expect_equal(e_xy(lt_bal, 60, 61, type = "complete", frac = "Balducci",
                    cohort = "last"), last_bal, tolerance = 1e-10)
})

test_that("e_xy preserves vector recycling, symmetry and tidy schema", {
  lt <- data.frame(x = 60:64, lx = c(100, 80, 50, 20, 0))
  attr(lt, "frac") <- "CF"
  for (type in c("curtate", "complete")) {
    for (cohort in c("first", "last")) {
      out <- e_xy(lt, c(60, 61), 62, t = c(1, 1.5), type = type,
                   cohort = cohort, tidy = TRUE)
      expect_identical(names(out), c("x", "y", "t", "type", "frac", "cohort", "ex"))
      expect_equal(out$frac, rep("CF", 2))
      expect_equal(out$ex, c(
        e_xy(lt, 62, 60, t = 1, type = type, cohort = cohort),
        e_xy(lt, 62, 61, t = 1.5, type = type, cohort = cohort)))
    }
  }
})
