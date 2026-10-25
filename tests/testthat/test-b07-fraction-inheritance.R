test_that("two-life survival inherits CML as constant force", {
  lt <- data.frame(x = 60:64, lx = 100 * exp(-.8 * (0:4)))
  attr(lt, "frac") <- "CML"
  other <- data.frame(x = 60:64, lx = 100 * exp(-.3 * (0:4)))
  attr(other, "frac") <- "CF"
  for (time in c(.25, .5, 1.5)) {
    sx <- exp(-.8 * time)
    sy <- exp(-.3 * time)
    expect_equal(t_pxy(lt, 60, 61, time), sx^2)
    expect_equal(t_pxy(lt, 60, 61, time, status = "last"), 2*sx - sx^2)
    expect_equal(t_pxy(list(lt, other), 60, 61, time), sx*sy)
    expect_equal(t_pxy(list(other, lt), 61, 60, time, status = "last"),
                 sx + sy - sx*sy)
    expect_equal(t_pxy(lt, 60, 61, time, frac = "UDD"),
                 (exp(-.8*floor(time)) *
                    (1-(time-floor(time))*(-expm1(-.8))))^2)
  }
  attr(other, "frac") <- "Balducci"
  expect_error(t_pxy(list(lt, other), 60, 60, .5), "different")
  expect_equal(t_pxy(list(lt, other), 60, 60, .5, frac = "CF"), exp(-.55))
})

test_that("two-life expectancy inherits CML with analytic exponential integrals", {
  mu <- .8
  lt <- data.frame(x = 60:64, lx = 100 * exp(-mu * (0:4)))
  attr(lt, "frac") <- "CML"
  for (horizon in c(.25, 1, 2.5)) {
    joint <- -expm1(-2*mu*horizon)/(2*mu)
    last <- 2*(-expm1(-mu*horizon)/mu) - joint
    expect_equal(e_xy(lt, 60, 61, t = horizon, type = "complete"), joint)
    expect_equal(e_xy(lt, 60, 61, t = horizon, type = "complete",
                      cohort = "last"), last)
  }
  out <- e_xy(lt, c(60, 61), 60, t = c(.5, 1), type = "complete", tidy = TRUE)
  expect_identical(out$frac, c("CF", "CF"))
  expect_equal(out$ex, -expm1(-2*mu*c(.5, 1))/(2*mu))
  expect_equal(e_xy(lt, 60, 61, t = 2), sum(exp(-2*mu*(1:2))))
})
