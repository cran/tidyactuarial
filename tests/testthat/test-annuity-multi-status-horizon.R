test_that("status matches legacy cohort and preserves positional calls", {
  lt <- data.frame(x = 60:64, lx = c(100, 80, 40, 10, 0))
  p <- c(.8, .4, .1)
  for (status in c("joint", "last")) {
    cohort <- if (status == "joint") "first" else "last"
    expected <- sum(if (status == "joint") p^2 else 2 * p - p^2)
    expect_equal(annuity_multi(lt, c(60, 60), 0, n = 3, status = status), expected)
    expect_equal(annuity_multi(lt, c(60, 60), 0, "effective", 1, 3, 0, 1,
                              "cohort", cohort, NULL, "immediate", "none", "UDD"),
                 expected)
  }
  expect_error(annuity_multi(lt, 60, 0, n = 2, cohort = "first", status = "joint"),
               "only one")
  expect_error(annuity_multi(lt, 60, 0, n = 2, status = "first"), "arg")
  expect_error(annuity_multi(lt, 60, 0, n = 2, annuity = "reversionary",
                             alpha = .5, status = "last"), "status.*reversionary")
})

test_that("common support uses years remaining from each life age", {
  a <- data.frame(x = 60:64, lx = c(100, 80, 40, 10, 0))
  b <- data.frame(x = 70:74, lx = c(200, 180, 120, 50, 0))
  tables <- list(a, b)
  expect_equal(annuity_multi(tables, c(60, 70), 0, n = 3), .985)
  expect_equal(annuity_multi(tables, c(60, 70), 0, n = 3, cohort = "last"), 2.065)
  expect_equal(annuity_multi(tables, c(60, 70), 0), .985)
  expect_equal(annuity_multi(rev(tables), c(70, 60), 0), .985)
  expect_equal(annuity_multi(tables, c(60, 70), 0, n = 2, h = 1,
                             annuity = "reversionary", cohort = "last", alpha = .5), .675)
  expect_error(annuity_multi(tables, c(60, 70), 0, n = 5), "horizon")
  expect_equal(annuity_multi(tables, c(60, 70), 0, n = 0, h = 5), 0)
  # The third table is open and short; even a last-survivor annuity must
  # stay within its observed support. At year one, P(any) = 1 - .2*.1*.3.
  c_table <- data.frame(x = 80:81, lx = c(100, 70))
  tables <- c(tables, list(c_table))
  expect_equal(annuity_multi(tables, c(60, 70, 80), 0, cohort = "last"), .994)
  expect_equal(annuity_multi(tables, c(60, 70, 80), 0, h = 1), 0)
  expect_error(annuity_multi(tables, c(60, 70, 80), 0, n = 2, cohort = "last"),
               "horizon")
})

test_that("unequal table ages retain the payment grid and Woolhouse borders", {
  a <- data.frame(x = 60:64, lx = 100 * exp(-.1 * (0:4)))
  b <- data.frame(x = 70:75, lx = 100 * exp(-.2 * (0:5)))
  tables <- list(a, b)
  F <- function(t) exp(-.35 * t)
  k <- 4
  for (timing in c("due", "immediate")) {
    times <- 1 + (if (timing == "due") 0:7 else 1:8) / k
    loss <- F(1) - F(3)
    first <- sum(F(1:2)) - (1 - 1 / k) * loss / 2
    if (timing == "immediate") first <- first - loss / k
    expected <- c(none = sum(F(times)) / k, first = first,
                  second = first - (1 - 1 / k^2) * .35 * loss / 12)
    for (order in names(expected)) {
      expect_equal(annuity_multi(tables, c(60, 70), .05, i_type = "force",
                                 n = 2, h = 1, k = k, timing = timing,
                                 frac = "CF", woolhouse = order, status = "joint"),
                   unname(expected[order]), tolerance = 1e-12)
    }
  }
  expect_error(annuity_multi(tables, c(60, 70), 0, n = 4, k = 2,
                             woolhouse = "second"), "Woolhouse")
})
