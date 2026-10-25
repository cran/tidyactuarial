test_that("annual mortality integrates each continuous force independently of frac", {
  cases <- list(
    list("Exponential", list(lambda = 0.1), function(s) rep(0.1, length(s))),
    list("Gompertz", list(B = 0.002, c = 1.09), function(s) 0.002 * 1.09^s),
    list("Makeham", list(A = 0.01, B = 0.002, c = 0.95),
         function(s) 0.01 + 0.002 * 0.95^s),
    list("Weibull", list(shape = 2.5, scale = 90),
         function(s) (2.5 / 90) * (s / 90)^1.5),
    list("Logistic", list(A = 0.01, B = 0.002, c = 1.09, C = 0.03),
         function(s) (0.01 + 0.002 * 1.09^s) / (1 + 0.03 * 1.09^s))
  )
  ages <- 40:43
  for (case in cases) {
    hazards <- vapply(ages, function(x) {
      integrate(case[[3]], x, x + 1, rel.tol = 1e-12)$value
    }, numeric(1))
    for (frac in c("CF", "CML", "UDD", "Balducci")) {
      tab <- mortality_law_table(case[[1]], 40, 43, params = case[[2]],
                                 frac = frac, close = FALSE, l0 = 1000)
      expect_equal(tab$qx, -expm1(-hazards), tolerance = 1e-12)
      expect_equal(tab$mu_x, case[[3]](ages), tolerance = 1e-12)
      expect_equal(tab$lx, 1000 * exp(-c(0, cumsum(hazards)[1:3])),
                   tolerance = 1e-12)
    }
  }
})

test_that("integrated forces handle constant-law limits and tiny probabilities", {
  for (c0 in c(1, 1 + 1e-12, 1 - 1e-12)) {
    tab <- mortality_law_table("Gompertz", 40, 42, B = 0.1, c = c0,
                               close = FALSE)
    ref <- vapply(40:42, function(x) {
      integrate(function(s) 0.1 * c0^s, x, x + 1)$value
    }, numeric(1))
    expect_equal(tab$qx, -expm1(-ref), tolerance = 1e-12)
  }
  tiny <- mortality_law_table("Exponential", 0, 2, lambda = 1e-20,
                              close = FALSE)
  expect_true(all(tiny$qx > 0))
  expect_equal(tiny$qx / 1e-20, rep(1, 3), tolerance = 1e-12)
  closed <- mortality_law_table("Gompertz", 40, 42, B = 0.01, c = 1.1)
  expect_equal(closed$qx[3], 1)
})
