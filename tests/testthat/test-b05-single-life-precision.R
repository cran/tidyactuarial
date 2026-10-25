test_that("single-life fractions retain observed small annual survival", {
  for (p in c(.8, .01, 1e-20, 1e-100, 1e-200, 1e-300)) {
    lt <- data.frame(x = 60:61, lx = c(1, p))
    for (s in c(.25, .5, .75)) {
      expect_equal(t_px(lt, 60, s, frac = "CF") / p^s, 1,
                   tolerance = 1e-12)
      expect_equal(t_px(lt, 60, s, frac = "CML") / p^s, 1,
                   tolerance = 1e-12)
      expect_equal(t_px(lt, 60, s, frac = "Balducci") /
                     (p / (p + (1-p)*s)), 1, tolerance = 1e-12)
    }
    for (s in c(.25, 1)) {
      mu <- -log(p)
      expect_equal(e_x(lt, 60, s, type = "complete", frac = "CF") /
                     (-expm1(-mu*s)/mu), 1, tolerance = 1e-12)
      expect_equal(e_x(lt, 60, s, type = "complete", frac = "Balducci") /
                     (p/(1-p)*log1p((1-p)*s/p)), 1, tolerance = 1e-12)
    }
  }
})

test_that("log survival resolves ratios outside the floating-point range", {
  lt <- data.frame(x = 60:61, lx = c(1e300, 1e-300))
  mu <- 600*log(10)
  expect_equal(t_px(lt, 60, .5, frac = "CF")/1e-300, 1, tolerance = 1e-12)
  expect_equal(e_x(lt, 60, 1, type = "complete", frac = "CF")*mu,
               1, tolerance = 1e-12)
  # Tiny duration compensates for an unrepresentable annual ratio.
  expect_equal(t_px(lt, 60, 1e-310, frac = "Balducci")/1e-290,
               1, tolerance = 1e-12)
  expect_equal(t_px(lt, 60, 1, frac = "CF"), 0)
})

test_that("precision fixes preserve columns, boundaries, NA and consumers", {
  lt <- data.frame(x = 60:63, lx = c(1, .5, 5e-201, 0))
  for (frac in c("UDD", "CF", "Balducci")) {
    expected <- switch(frac, UDD = .75, CF = sqrt(.5), Balducci = 2/3)
    expect_equal(t_px(lt, 60, .5, frac = frac), expected)
    expect_equal(t_px(lt, 60, c(0, 1, 3, 4), frac = frac), c(1,.5,0,0))
    expect_equal(e_x(lt, 60, t = 1, frac = frac), .5)
    # Explicit closed forms for each interval, independent of t_px/e_x.
    area <- function(p, s) switch(frac, UDD = s - (1-p)*s*s/2,
      CF = -expm1(log(p)*s)/(-log(p)),
      Balducci = p/(1-p)*log1p((1-p)*s/p))
    expect_equal(e_x(lt, 60, 1.5, type = "complete", frac = frac),
                 area(.5,1) + .5*area(1e-200,.5), tolerance = 1e-12)
    expect_equal(t_qx(lt, 60, .5, frac = frac), 1-expected)
  }
  tiny <- data.frame(x = 60:61, lx = c(1, 1e-200))
  with_px <- transform(tiny, px = c(1e-200, NA_real_))
  with_qx <- transform(tiny, qx = c(1, NA_real_))
  both <- transform(with_px, qx = c(1, NA_real_))
  for (tab in list(tiny, with_px, with_qx, both)) {
    expect_equal(t_px(tab, 60, .5, frac = "CF")/1e-100, 1)
    expect_equal(e_x(tab, 60, 1, type = "complete", frac = "CF") *
                   (200*log(10)), 1)
  }
  # Supplied px remains authoritative even when it differs from lx.
  explicit <- transform(tiny, px = c(.25, NA_real_), qx = c(.75, NA_real_))
  expect_equal(t_px(explicit, 60, .5, frac = "CF"), .5)
  expect_equal(e_x(explicit, 60, .5, type = "complete", frac = "CF"),
               .5/log(4))
  expect_named(t_px(tiny, 60, .5, tidy = TRUE), c("x","t","frac","tpx"))
  expect_named(e_x(tiny, 60, .5, tidy = TRUE), c("x","t","type","frac","ex"))
  missing <- data.frame(x = 60:62, lx = c(1,NA_real_,0))
  expect_true(is.na(t_px(missing, 60, .5, frac = "CF")))
  expect_true(is.na(e_x(missing, 60, type = "complete", frac = "CF")))
  extinct <- data.frame(x = 60:62, lx = c(1,0,0))
  expect_equal(e_x(extinct, 60, type = "complete", frac = "CF"), 0)
  expect_equal(e_x(extinct, 60, type = "complete", frac = "UDD"), .5)
})
