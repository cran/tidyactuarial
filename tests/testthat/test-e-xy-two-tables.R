test_that("two-table curtate expectations equal independent survival ledgers", {
  a <- data.frame(x = 60:63, lx = c(100, 80, 40, 10))
  b <- data.frame(x = 70:75, lx = c(200, 180, 140, 80, 40, 20))
  px <- c(.8, .4, .1, 0, 0)
  py <- c(.9, .7, .4, .2, .1)
  for (n in c(0, .5, 1, 2.5, 3, 5, 9)) {
    keep <- seq_along(px) <= floor(n)
    expect_equal(e_xy(list(a, b), 60, 70, t = n, status = "joint"),
                 sum((px * py)[keep]))
    expect_equal(e_xy(list(a, b), 60, 70, t = n, status = "last"),
                 sum((px + py - px * py)[keep]))
  }
  expect_equal(e_xy(list(a, b), 60, 70), sum(px * py))
  expect_equal(e_xy(list(a, b), 60, 70, status = "last"),
               sum(px + py - px * py))
  # The originally reported closed-table example has only one nonzero term.
  a <- data.frame(x = 60:62, lx = c(100, 80, 0))
  b <- data.frame(x = 60:62, lx = c(100, 50, 0))
  expect_equal(e_xy(list(a, b), 60, 60, status = "joint"), .4)
  expect_equal(e_xy(list(a, b), 60, 60, status = "last"), .9)
})

test_that("complete two-table UDD expectations integrate the original marginals", {
  a <- data.frame(x = 60:63, lx = c(100, 80, 40, 10))
  b <- data.frame(x = 70:75, lx = c(200, 180, 140, 80, 40, 20))
  # Independent piecewise-linear marginals, zero beyond each observed horizon.
  sx <- function(u) ifelse(u > 3, 0, approx(0:3, a$lx / 100, u, rule = 2)$y)
  sy <- function(u) ifelse(u > 5, 0, approx(0:5, b$lx / 200, u, rule = 2)$y)
  for (n in c(0, .3, 1, 2.5, 4.5, 5, 8)) {
    end <- min(n, 5)
    cuts <- sort(unique(c(0, (0:5)[0:5 < end], end)))
    joint <- last <- 0
    if (end > 0) for (j in seq_len(length(cuts) - 1L)) {
      joint <- joint + integrate(function(u) sx(u) * sy(u),
                                  cuts[j], cuts[j + 1L])$value
      last <- last + integrate(function(u) sx(u) + sy(u) - sx(u) * sy(u),
                                cuts[j], cuts[j + 1L])$value
    }
    for (status in c("joint", "last")) {
      got <- e_xy(list(a, b), 60, 70, t = n, type = "complete",
                  frac = "UDD", status = status)
      expect_equal(got, if (status == "joint") joint else last, tolerance = 1e-12)
    }
  }
})

test_that("unequal-table CF and Balducci expectations match analytic integrals", {
  for (frac in c("CF", "CML", "Balducci")) {
    for (n in c(.7, 2.5, 8)) {
      hx <- min(n, 3); hy <- min(n, 5); hj <- min(hx, hy)
      if (frac %in% c("CF", "CML")) {
        a <- data.frame(x = 60:63, lx = 100 * exp(-.2 * (0:3)))
        b <- data.frame(x = 70:75, lx = 200 * exp(-.35 * (0:5)))
        joint <- -expm1(-.55 * hj) / .55
        total <- -expm1(-.2 * hx) / .2 - expm1(-.35 * hy) / .35
      } else {
        a <- data.frame(x = 60:63, lx = 100 * 2 / (2 + 0:3))
        b <- data.frame(x = 70:75, lx = 200 * 5 / (5 + 0:5))
        joint <- 10 / 3 * (log1p(hj / 2) - log1p(hj / 5))
        total <- 2 * log1p(hx / 2) + 5 * log1p(hy / 5)
      }
      got_joint <- e_xy(list(a, b), 60, 70, t = n, type = "complete",
                        frac = frac, status = "joint")
      got_last <- e_xy(list(a, b), 60, 70, t = n, type = "complete",
                       frac = frac, status = "last")
      expect_equal(got_joint, joint, tolerance = 1e-12)
      expect_equal(got_last, total - joint, tolerance = 1e-12)
      expect_equal(got_joint + got_last,
                   e_x(a, 60, t = n, type = "complete", frac = frac) +
                     e_x(b, 70, t = n, type = "complete", frac = frac),
                   tolerance = 1e-12)
    }
  }
})

test_that("status preserves legacy calls and the existing tidy schema", {
  a <- data.frame(x = 60:64, lx = c(100, 80, 50, 20, 0))
  for (type in c("curtate", "complete")) {
    for (status in c("joint", "last")) {
      cohort <- if (status == "joint") "first" else "last"
      old <- e_xy(a, c(60, 61), 62, c(1, 1.5), type, "CF", cohort, TRUE, TRUE, 1e-10)
      new <- e_xy(a, c(60, 61), 62, c(1, 1.5), type, "CF", tidy = TRUE, status = status)
      expect_identical(new, old)
      expect_identical(names(new), c("x", "y", "t", "type", "frac", "cohort", "ex"))
      expect_identical(e_xy(list(a, a), 60, 61, type = type, status = status),
                       e_xy(a, 60, 61, type = type, cohort = cohort))
    }
  }
  expect_error(e_xy(a, 60, 61, cohort = "first", status = "joint"),
               "only one of.*status.*cohort")
  expect_error(e_xy(a, 60, 61, status = "first"), "arg")
  expect_error(e_xy(a, 60, 61, status = NA_character_), "arg")
})

test_that("two-table fractional assumptions follow t_pxy inheritance rules", {
  a <- data.frame(x = 60:62, lx = c(100, 80, 20))
  b <- data.frame(x = 70:73, lx = c(100, 70, 30, 10))
  attr(a, "frac") <- "CML"
  for (second in c("CF", "CML")) {
    attr(b, "frac") <- second
    expect_equal(e_xy(list(a, b), 60, 70, type = "complete"),
                 e_xy(list(a, b), 60, 70, type = "complete", frac = "CF"))
  }
  attr(b, "frac") <- NULL
  expect_equal(e_xy(list(b, a), 70, 60, type = "complete"),
               e_xy(list(a, b), 60, 70, type = "complete", frac = "CF"))
  attr(b, "frac") <- "Balducci"
  expect_error(e_xy(list(a, b), 60, 70), "different.*frac.*explicitly")
  expect_equal(e_xy(list(a, b), 60, 70, type = "complete", frac = "UDD"),
               e_xy(list(b, a), 70, 60, type = "complete", frac = "UDD"))
})

test_that("each table is validated and each initial life uses its own marginal", {
  a <- data.frame(x = 60:63, lx = c(100, 80, 40, 10))
  b <- data.frame(x = 70:75, lx = c(200, 180, 140, 80, 40, 0))
  for (status in c("joint", "last")) for (type in c("curtate", "complete")) {
    got <- e_xy(list(a, b), c(60, 61), 71, t = c(1, 2.5),
                type = type, status = status)
    swapped <- e_xy(list(b, a), 71, c(60, 61), t = c(1, 2.5),
                    type = type, status = status)
    expect_equal(got, swapped)
    expect_true(is.na(e_xy(list(a, b), 60, 75, t = 0, type = type, status = status)))
    expect_true(is.na(e_xy(list(a, b), 70, 60, t = 0, type = type, status = status)))
    expect_equal(e_xy(list(a, b), 60, 70, t = 0, type = type, status = status), 0)
  }
  # Positive initial lx at its final recorded age: zero complete joint area;
  # last-survivor expectation still contains the other marginal's future.
  for (type in c("curtate", "complete")) {
    expect_equal(e_xy(list(a, b), 63, 70, type = type, status = "joint"), 0)
    expect_equal(e_xy(list(a, b), 63, 70, type = type, status = "last"),
                 e_x(b, 70, type = type))
  }
  expect_error(e_xy(list(a, b[-2, ]), 60, 70, check = FALSE), "consecutive")
  expect_error(e_xy(list(a, transform(b, lx = rev(lx))), 60, 70), "nonincreasing")
  expect_error(e_xy(list(a), 60, 70), "one life table.*two life tables")
  expect_error(e_xy(list(a, 1), 60, 70), "one life table.*two life tables")
  expect_error(e_xy(list(a, b["x"]), 60, 70), "columns")
})
