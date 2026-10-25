# Independent exposure calculations: Lx = lx_next + ax * dx, mx = dx / Lx.
# Lx is a test reference, not an additional public output column.
for (a in c(0, 0.5, 1)) {
  local({
    ax_value <- a
    test_that(paste("lifetable respects exposure for ax =", ax_value), {
      survivors <- c(1000, 800, 400)
      deaths <- c(200, 400, 400)
      exposure <- c(800, 400, 0) + ax_value * deaths
      reference_m <- deaths / exposure
      from_l <- lifetable(60:62, lx = survivors, ax = ax_value)
      from_q <- lifetable(60:62, qx = c(0.2, 0.5, 1),
                          radix = 1000, ax = ax_value)
      for (lt in list(from_l, from_q)) {
        expect_equal(lt$lx, survivors)
        expect_equal(lt$dx, deaths)
        expect_equal(lt$qx, deaths / survivors)
        expect_equal(lt$px, c(0.8, 0.5, 0))
        expect_equal(lt$mx, reference_m, tolerance = 1e-12)
        expect_equal(lt$mx[1:2] * exposure[1:2], deaths[1:2])
        expect_identical(names(lt), c("x", "lx", "dx", "qx", "px", "mx"))
      }
      # Keep the documented rejection of infinite mx inputs. The infinite
      # terminal output at ax = 0 is tested through lx and qx above.
      inverse <- lifetable(60:62, mx = c(reference_m[1:2], NA_real_),
                           radix = 1000, ax = ax_value)
      expect_equal(inverse$lx, survivors, tolerance = 1e-12)
      expect_equal(inverse$dx, deaths, tolerance = 1e-12)
      expect_equal(inverse$qx, deaths / survivors, tolerance = 1e-12)
      expect_equal(inverse$mx, reference_m, tolerance = 1e-12)
    })
  })
}

test_that("constant force agrees with independently integrated exposure", {
  mu <- 0.1
  q <- -expm1(-mu)
  # Conditional expected death time and person-years from continuous survival.
  a <- integrate(function(t) t * mu * exp(-mu * t), 0, 1)$value / q
  ages <- 60:63
  survivors <- 1000 * exp(-mu * (0:3))
  exposure <- vapply(survivors, function(l) {
    integrate(function(t) l * exp(-mu * t), 0, 1)$value
  }, numeric(1))
  deaths <- survivors * q
  from_q <- lifetable(ages, qx = rep(q, 4), radix = 1000,
                      close = FALSE, ax = a, frac = "CF")
  from_m <- lifetable(ages, mx = rep(mu, 4), radix = 1000,
                      close = FALSE, ax = a, frac = "CF")
  expect_equal(a, 1 / mu - 1 / expm1(mu), tolerance = 1e-12)
  expect_equal(exposure, survivors - deaths + a * deaths, tolerance = 1e-12)
  for (lt in list(from_q, from_m)) {
    expect_equal(lt$lx, survivors, tolerance = 1e-12)
    expect_equal(lt$dx, deaths, tolerance = 1e-12)
    expect_equal(lt$qx, rep(q, 4), tolerance = 1e-12)
    expect_equal(lt$mx, deaths / exposure, tolerance = 1e-12)
    expect_equal(lt$mx, rep(mu, 4), tolerance = 1e-12)
  }
})

test_that("terminal mx uses the exact ax limit for positive survivors", {
  for (a in c(0, 1e-12, 0.25, 0.5, 0.75, 1)) {
    for (closed in c(FALSE, TRUE)) {
      lt <- lifetable(60:61, qx = c(0.2, 1), radix = 1000,
                      ax = a, close = closed)
      expect_equal(lt$lx[2], 800)
      expect_equal(lt$dx[2], 800)
      expect_equal(lt$qx[2], 1)
      # L_omega = a * 800. At a = 0 the positive-deaths/zero-exposure
      # ratio has the extended-real limit +Inf; at a = 1 it equals 1.
      expect_equal(lt$mx[2], 800 / (a * 800), tolerance = 1e-12)
    }
  }
  for (a in c(0.25, 0.5, 0.75, 1)) {
    lt <- lifetable(60:61, mx = c(0.1, 1 / a), radix = 1000,
                    ax = a, close = FALSE)
    expect_equal(lt$qx[2], 1, tolerance = 1e-12)
    expect_equal(lt$dx[2], lt$lx[2], tolerance = 1e-12)
  }
  expect_error(lifetable(60:61, mx = c(0.1, Inf), radix = 1000, ax = 0),
               "must not contain Inf")
})

test_that("UDD results preserve the existing ax = 0.5 contract", {
  lt <- lifetable(60:62, lx = c(1000, 800, 400))
  expect_equal(lt$mx, c(2 / 9, 2 / 3, 2), tolerance = 1e-12)
  expect_equal(lt, lifetable(60:62, lx = c(1000, 800, 400), ax = 0.5))
  inverse <- lifetable(60:62, mx = c(2 / 9, 2 / 3, 2), radix = 1000)
  expect_equal(inverse, lt, tolerance = 1e-12)
})

test_that("actuarial consumers agree for equivalent qx and mx tables", {
  survivors <- c(1000, 800, 400, 200)
  deaths <- c(200, 400, 200, 200)
  a <- 0.25
  exposure <- c(800, 400, 200, 0) + a * deaths
  from_q <- lifetable(60:63, qx = deaths / survivors, radix = 1000, ax = a)
  from_m <- lifetable(60:63, mx = c((deaths / exposure)[1:3], NA_real_),
                      radix = 1000, ax = a)
  expect_equal(t_px(from_m, x = 60, t = 2), 0.4, tolerance = 1e-12)
  expect_equal(e_x(from_m, x = 60), e_x(from_q, x = 60), tolerance = 1e-12)
  expect_equal(annuity_x(from_m, x = 60, i = 0.03),
               annuity_x(from_q, x = 60, i = 0.03), tolerance = 1e-12)
  expect_equal(insurance_x(from_m, x = 60, i = 0.03),
               insurance_x(from_q, x = 60, i = 0.03), tolerance = 1e-12)
  expect_equal(premium_x(from_m, x = 60, i = 0.03),
               premium_x(from_q, x = 60, i = 0.03), tolerance = 1e-12)
})
