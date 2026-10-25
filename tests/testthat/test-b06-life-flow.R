test_that("life flows discount from the original interest specification", {
  lt <- data.frame(x = 60:62, lx = c(100, 100, 100))
  for (delta in c(-1000, 1000)) {
    out <- apv_life_flow(lt, 60, t = .001, cf = 2,
                         i = delta, i_type = "force")
    expect_equal(out$discount, exp(-delta * .001))
    expect_equal(out$pv, 2 * exp(-delta * .001))
  }
  out <- apv_life_flow(lt, 60, t = 1, cf = 1, i = 1e-17)
  expect_equal(out$discount, exp(-log1p(1e-17)))
  out <- apv_life_flow(lt, 60, t = .01, cf = 2,
                       i = 100, i_type = "nominal_interest", m = 100)
  expect_equal(out$pv, 1) # 100 conversion periods per year, one period here.
  out <- apv_life_flow(lt, 60, t = .01, cf = 2,
                       i = 50, i_type = "nominal_discount", m = 100)
  expect_equal(out$pv, 1)
})

test_that("life-flow monetary products retain range and signs", {
  lt <- data.frame(x = 60:61, px = c(.5, .5))
  out <- apv_life_flow(lt, 60, t = c(1, 1, 1), cf = c(1e300, -1e300, 0),
                       i = 800, i_type = "force")
  reference <- exp(log(1e300) + log(.5) - 800)
  expect_equal(out$pv / reference, c(1, -1, 0))
  expect_identical(out$discount, rep(0, 3))
  out <- apv_life_flow(lt, 60, t = 1, cf = 1e-300,
                       i = -800, i_type = "force")
  expect_equal(out$pv / exp(log(1e-300) + log(.5) + 800), 1)
  expect_identical(out$discount, Inf)
  lt <- data.frame(x = 60:61, px = c(1e-200, 1))
  out <- apv_life_flow(lt, 60, t = 1, cf = 1e-200,
                       i = -700, i_type = "force")
  expect_equal(out$pv / exp(2 * log(1e-200) + 700), 1)
  expect_identical(out$expected_cf, 0)
})

test_that("missing life-flow times propagate to rows without losing positions", {
  lt <- data.frame(x = 60:62, lx = c(100, 80, 0))
  out <- apv_life_flow(lt, 60, t = c(0, NA, 1), cf = c(10, 20, 30), i = 0)
  expect_equal(out$surv_prob, c(1, NA, .8))
  expect_equal(out$pv, c(10, NA, 24))
  expect_equal(out$pv_cum, c(10, NA, NA))
  expect_equal(attr(out, "apv"), 34) # Existing nonmissing-subtotal convention.
  out <- apv_life_flow(lt, 60, t = 0:1, cf = c(NA, 10), i = 0)
  expect_equal(out$pv, c(NA, 8))
  expect_equal(attr(out, "apv"), 8)
})

test_that("last-survivor flow preserves rare independent survival", {
  lt <- data.frame(x = 60:62, px = c(1e-20, 2e-20, .5))
  out <- apv_life_flow(lt, c(60, 61), t = 1, cf = 1,
                       i = 0, status = "last")
  expect_equal(out$surv_prob / 3e-20, 1)
  out <- apv_life_flow(lt, c(60, 61), t = 1, cf = 1,
                       i = 0, status = "reversionary", alpha = .4)
  expect_equal(out$surv_prob / 1.2e-20, 1)
})

test_that("life-flow horizon queries are bounded by available mortality rows", {
  lt <- data.frame(x = 60:62, px = c(.8, 0, NA_real_))
  out <- apv_life_flow(lt, 60, t = c(1, 2, 2.5, 1e10), cf = rep(1, 4), i = 0)
  expect_equal(out$surv_prob, c(.8, 0, 0, 0))
  lt$px[2] <- .5
  out <- apv_life_flow(lt, 60, t = c(2, 3.5, 1e10), cf = rep(1, 3), i = 0)
  expect_equal(out$surv_prob, c(.4, NA, NA))
})

test_that("life-flow integer and date inputs avoid coercion failures", {
  lt <- data.frame(x = 60:62, px = c(.8, .5, 0))
  expect_error(apv_life_flow(lt, 1e20, t = 0, cf = 1, i = 0), "ages")
  expect_error(apv_life_flow(lt, 60, date = as.Date(Inf, origin = "1970-01-01"),
                             cf = 1, i = 0), "date")
  expect_error(apv_life_flow(lt, 60, date = as.Date("2026-01-01"),
                             date0 = as.Date(-Inf, origin = "1970-01-01"),
                             cf = 1, i = 0), "date0")
  expect_no_warning(out <- apv_life_flow(lt, 60, t = 1, cf = 1, i = .1, m = 3e9))
  expect_equal(out$pv, .8 / 1.1)
})

test_that("life-flow mortality validation prevents invalid probabilities", {
  bad <- list(
    data.frame(x = c(60, 60), px = c(.8, .5)),
    data.frame(x = c(60, Inf), px = c(.8, .5)),
    data.frame(x = c(60, 61.5), px = c(.8, .5)),
    data.frame(x = 60:61, px = c(1.1, .5)),
    data.frame(x = 60:61, qx = c(-.1, .5)),
    data.frame(x = 60:61, lx = c(100, 110)),
    data.frame(x = 60:61, lx = c(100, Inf)),
    data.frame(x = 60:61, px = c(.8, NaN)),
    data.frame(x = 60:61, px = c(.8, .5), qx = c(.3, .5))
  )
  for (lt in bad) {
    expect_error(apv_life_flow(lt, 60, t = .5, cf = 1, i = 0), "lt")
  }
  # Gaps/missing mortality are unknown, not silently interpolated or dropped.
  lt <- data.frame(x = c(60, 62), px = c(.8, .5))
  out <- apv_life_flow(lt, 60, t = 2, cf = 1, i = 0)
  expect_true(is.na(out$pv))
  out <- apv_life_flow(lt, 60, t = numeric(), cf = numeric(), i = 0)
  expect_equal(nrow(out), 0)
  expect_identical(attr(out, "apv"), 0)
})

test_that("life-flow UDD and multi-life formulas retain independent references", {
  lt <- data.frame(x = 60:63, px = c(.8, .5, .25, 0))
  t <- c(1.5, 0, 1.5, .5)
  cf <- c(-100, 10, 200, 30)
  single <- c(.8 * .75, 1, .8 * .75, .9)
  second <- c(.5 * .625, 1, .5 * .625, .75)
  delta <- log1p(.05)
  for (status in c("single", "first", "last", "reversionary")) {
    p <- switch(status, single = single, first = single * second,
                last = single + second - single * second,
                reversionary = single * second + .4 *
                  (single + second - 2 * single * second))
    out <- apv_life_flow(lt, c(60, 61), t = t, cf = cf,
                         i = .05, status = status, alpha = .4)
    expect_equal(out$t, t)
    expect_equal(out$surv_prob, p)
    expect_equal(out$pv, cf * p * exp(-delta * t))
    expect_equal(out$pv_cum, cumsum(cf * p * exp(-delta * t)))
    expect_equal(attr(out, "apv"), sum(cf * p * exp(-delta * t)))
  }
})
