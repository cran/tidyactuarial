test_that("IRR candidates are invariant to monetary units", {
  for (scale in c(1, 1e-9, 1e-200, 1e200)) {
    out <- irr_flow_multi(cf = scale * c(-100, 110), t = c(0, 1),
      search_interval = c(0, .5), grid_points = 101)
    expect_equal(nrow(out), 1L)
    expect_equal(out$irr, .1, tolerance = 1e-9)
    one <- irr_flow(cf = scale * c(-100, 110), t = c(0, 1),
      interval = c(0, .5))
    expect_true(one$converged)
    expect_equal(one$irr, .1, tolerance = 1e-9)
    # Residual stays in the original currency, not in normalized units.
    expected <- sum(c(-100, 110) / (1 + one$irr)^c(0, 1))
    expect_lt(abs(one$npv / scale - expected), 1e-10)
  }
})

test_that("tiny cash flows do not make every grid point an IRR", {
  # -100 + 50/(1+r) is strictly negative throughout [0,1].
  out <- irr_flow_multi(cf = 1e-200 * c(-100, 50), t = c(0, 1),
    search_interval = c(0, 1), grid_points = 101)
  expect_equal(nrow(out), 0L)
  one <- irr_flow(cf = 1e-200 * c(-100, 50), t = c(0, 1),
    interval = c(0, 1))
  expect_false(one$converged)
})

test_that("scaled quadratic cash flows retain two independently known roots", {
  # 50(1+r)^2 - 115(1+r) + 66 = 50(r-.1)(r-.2).
  for (scale in c(1e-9, 1e-200, 1e200)) {
    out <- irr_flow_multi(cf = scale * c(50, -115, 66), t = 0:2,
      search_interval = c(0, .5), grid_points = 101)
    expect_equal(nrow(out), 2L)
    expect_equal(out$irr, c(.1, .2), tolerance = 1e-8)
  }
})

test_that("IRR scaling rejects a loss of nonzero input amounts", {
  for (fun in list(irr_flow, irr_flow_multi)) {
    expect_error(fun(cf = c(1e308, -1e-308), t = c(0, 1)),
      "lose nonzero amounts")
  }
})
