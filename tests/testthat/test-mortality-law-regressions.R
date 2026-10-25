test_that("Heligman-Pollard uses eight parameters with separate components", {
  # At x = F the accident component is exactly D. Each component can
  # also be isolated by setting the other amplitudes to zero.
  pars <- list(A = 0.01, B = 1, C = 0.5, D = 0.02, E = 2,
               F_hp = 3, G = 0.001, H = 2)
  tab <- mortality_law_table("HeligmanPollard", 3, 4, params = pars,
                             close = FALSE)
  odds_at_3 <- 0.01^2 + 0.02 + 0.001 * 8
  expect_equal(tab$qx[1], odds_at_3 / (1 + odds_at_3), tolerance = 1e-12)
  for (component in c("child", "accident", "old")) {
    isolated <- pars
    isolated$A <- if (component == "child") pars$A else 0
    isolated$D <- if (component == "accident") pars$D else 0
    isolated$G <- if (component == "old") pars$G else 0
    odds <- switch(component, child = 0.01^2, accident = 0.02, old = 0.008)
    result <- mortality_law_table("HeligmanPollard", 3, 4,
                                  params = isolated, close = FALSE)
    expect_equal(result$qx[1], odds / (1 + odds), tolerance = 1e-12)
  }
  missing_h <- pars[names(pars) != "H"]
  expect_error(mortality_law_table("HeligmanPollard", 3, 4,
                                   params = missing_h), "Missing parameter.*H")
  for (bad in list(list(F_hp = 0), list(H = 0), list(B = -4), list(C = 0))) {
    expect_error(mortality_law_table("HeligmanPollard", 3, 4,
                                     params = utils::modifyList(pars, bad)))
  }
})

test_that("Beta conditional mortality remains accurate in the upper tail", {
  # Beta(1, b): S(x) = (1 - x / omega)^b, even far below machine epsilon.
  tab <- mortality_law_table("Beta", 90, 99, alpha = 1, beta = 20,
                             omega = 100, close = FALSE)
  expected <- 1 - ((100 - (90:99) - 1) / (100 - (90:99)))^20
  expect_equal(tab$qx, expected, tolerance = 1e-12)
})

test_that("Heligman-Pollard remains finite when old-age odds overflow", {
  tab <- mortality_law_table("HeligmanPollard", 100, 102,
    A = 0, B = 0, C = 1, D = 0, E = 1, F_hp = 20,
    G = 1, H = 1e10, close = FALSE)
  expect_equal(tab$qx, rep(1, 3))
  expect_false(anyNA(tab$qx))
})

test_that("Heligman-Pollard rejects invalid parameters and ambiguous aliases", {
  pars <- list(A = 0.01, B = 1, C = 0.5, D = 0.02, E = 2,
               F_hp = 3, G = 0.001, H = 2)
  for (nm in names(pars)) {
    for (bad in c(NA_real_, Inf, -1)) {
      invalid <- pars
      invalid[[nm]] <- bad
      expect_error(mortality_law_table("HeligmanPollard", 3, 4, params = invalid))
    }
  }
  expect_error(mortality_law_table("HeligmanPollard", 3, 4,
    params = c(pars, list(F = 3))), "Provide only one")
  closed <- mortality_law_table("HeligmanPollard", 3, 4, params = pars)
  expect_equal(closed$qx[2], 1)
  expect_identical(names(closed), c("x", "law", "frac", "mu_x", "qx", "px",
                                   "lx", "dx", "Lx", "Tx", "ex", "mx"))
})

test_that("exposure ratios do not depend on the radix magnitude", {
  for (radix in c(1000, 1e-20)) {
    tab <- mortality_law_table("Exponential", 0, 2, lambda = 0.1,
                               a_x = 0, l0 = radix)
    expect_equal(tab$mx[1:2], rep(expm1(0.1), 2), tolerance = 1e-12)
    expect_identical(tab$mx[3], Inf)
    expect_equal(tab$ex, tab$Tx / tab$lx, tolerance = 1e-12)
  }
})

test_that("extinct cohorts retain undefined exposure ratios", {
  tab <- mortality_law_table("DeMoivre", 98, 100, omega = 100)
  expect_equal(tab$lx[3], 0)
  expect_true(is.nan(tab$mx[3]))
  expect_true(is.nan(tab$ex[3]))
})
