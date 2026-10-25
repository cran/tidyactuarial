test_that("mortality_law_table Exponential reproduces constant-force survival", {
  lambda <- 0.02

  tab <- mortality_law_table(
    law = "Exponential",
    x_min = 40,
    x_max = 43,
    lambda = lambda,
    frac = "CF",
    l0 = 1000,
    close = FALSE
  )

  q <- 1 - exp(-lambda)
  p <- exp(-lambda)

  expect_s3_class(tab, "tbl_df")
  expect_equal(tab$x, 40:43)
  expect_equal(tab$mu_x, rep(lambda, 4), tolerance = 1e-12)
  expect_equal(tab$qx, rep(q, 4), tolerance = 1e-12)
  expect_equal(tab$px, rep(p, 4), tolerance = 1e-12)
  expect_equal(
    tab$lx,
    1000 * p^(0:3),
    tolerance = 1e-10
  )
})

test_that("annual probabilities are independent of the interpolation label", {
  lambda <- 0.10

  cf <- mortality_law_table(
    "Exponential", 30, 32,
    lambda = lambda,
    frac = "CF",
    close = FALSE
  )

  udd <- mortality_law_table(
    "Exponential", 30, 32,
    lambda = lambda,
    frac = "UDD",
    close = FALSE
  )

  bal <- mortality_law_table(
    "Exponential", 30, 32,
    lambda = lambda,
    frac = "Balducci",
    close = FALSE
  )

  expect_equal(cf$qx, rep(1 - exp(-lambda), 3), tolerance = 1e-12)
  expect_equal(udd$qx, rep(1 - exp(-lambda), 3), tolerance = 1e-12)
  expect_equal(
    bal$qx,
    rep(1 - exp(-lambda), 3),
    tolerance = 1e-12
  )
})

test_that("CML is the documented alias for constant force", {
  cf <- mortality_law_table(
    "Gompertz", 50, 53,
    B = 1e-5,
    c = 1.08,
    frac = "CF",
    close = FALSE
  )

  cml <- mortality_law_table(
    "Gompertz", 50, 53,
    B = 1e-5,
    c = 1.08,
    frac = "CML",
    close = FALSE
  )

  expect_equal(cml$qx, cf$qx, tolerance = 1e-12)
  expect_true(all(cml$frac == "CF"))
})

test_that("mortality_law_table Gompertz and Makeham reproduce their forces", {
  ages <- 40:43

  B <- 2e-5
  c0 <- 1.09

  gompertz <- mortality_law_table(
    "Gompertz",
    min(ages),
    max(ages),
    B = B,
    c = c0,
    frac = "CF",
    close = FALSE
  )

  mu_g <- B * c0^ages

  expect_equal(gompertz$mu_x, mu_g, tolerance = 1e-12)
  expect_equal(gompertz$qx, -expm1(-mu_g * (c0 - 1) / log(c0)), tolerance = 1e-12)

  A <- 4e-4

  makeham <- mortality_law_table(
    "Makeham",
    min(ages),
    max(ages),
    A = A,
    B = B,
    c = c0,
    frac = "CF",
    close = FALSE
  )

  mu_m <- A + B * c0^ages

  expect_equal(makeham$mu_x, mu_m, tolerance = 1e-12)
  hazard_m <- A + (mu_m - A) * (c0 - 1) / log(c0)
  expect_equal(makeham$qx, -expm1(-hazard_m), tolerance = 1e-12)
})

test_that("mortality_law_table Weibull and Logistic reproduce their forces", {
  ages <- 20:23

  shape <- 2.5
  scale <- 90

  weibull <- mortality_law_table(
    "Weibull",
    min(ages),
    max(ages),
    shape = shape,
    scale = scale,
    frac = "CF",
    close = FALSE
  )

  mu_w <- (shape / scale) * ((ages / scale)^(shape - 1))

  expect_equal(weibull$mu_x, mu_w, tolerance = 1e-12)
  hazard_w <- ((ages + 1) / scale)^shape - (ages / scale)^shape
  expect_equal(weibull$qx, -expm1(-hazard_w), tolerance = 1e-12)

  A <- 1e-4
  B <- 2e-6
  c0 <- 1.10
  C <- 1e-3

  logistic <- mortality_law_table(
    "Logistic",
    min(ages),
    max(ages),
    A = A,
    B = B,
    c = c0,
    C = C,
    frac = "CF",
    close = FALSE
  )

  cx <- c0^ages
  mu_l <- (A + B * cx) / (1 + C * cx)

  expect_equal(logistic$mu_x, mu_l, tolerance = 1e-12)
  hazard_l <- vapply(ages, function(x) {
    integrate(function(s) (A + B * c0^s) / (1 + C * c0^s), x, x + 1)$value
  }, numeric(1))
  expect_equal(logistic$qx, -expm1(-hazard_l), tolerance = 1e-12)
})

test_that("Weibull transitional aliases agree with the current parameter names", {
  current <- mortality_law_table(
    "Weibull", 20, 24,
    shape = 2.5,
    scale = 90,
    close = FALSE
  )

  legacy <- mortality_law_table(
    "Weibull", 20, 24,
    k = 2.5,
    lambda = 90,
    close = FALSE
  )

  expect_equal(legacy$mu_x, current$mu_x, tolerance = 1e-12)
  expect_equal(legacy$qx, current$qx, tolerance = 1e-12)
})

test_that("De Moivre life table follows q_x = 1 / (omega - x)", {
  omega <- 65
  ages <- 60:63

  tab <- mortality_law_table(
    "DeMoivre",
    min(ages),
    max(ages),
    omega = omega,
    l0 = 1000,
    close = FALSE
  )

  expected_q <- 1 / (omega - ages)

  expect_true(all(is.na(tab$mu_x)))
  expect_true(all(is.na(tab$frac)))
  expect_equal(tab$qx, expected_q, tolerance = 1e-12)

  # Under De Moivre, successive survivor counts are proportional to omega - x.
  expected_lx <- 1000 * (omega - ages) / (omega - ages[1])
  expect_equal(tab$lx, expected_lx, tolerance = 1e-10)
})

test_that("Beta mortality law matches survival ratios from the beta distribution", {
  alpha <- 2
  beta <- 5
  omega <- 100
  ages <- 60:63

  tab <- mortality_law_table(
    "Beta",
    min(ages),
    max(ages),
    alpha = alpha,
    beta = beta,
    omega = omega,
    close = FALSE
  )

  Sx <- 1 - stats::pbeta(ages / omega, shape1 = alpha, shape2 = beta)
  Sx1 <- 1 - stats::pbeta((ages + 1) / omega, shape1 = alpha, shape2 = beta)

  expected_q <- (Sx - Sx1) / Sx

  expect_true(all(is.na(tab$mu_x)))
  expect_equal(tab$qx, expected_q, tolerance = 1e-12)
})

test_that("Heligman-Pollard mortality law reproduces its odds specification", {
  ages <- 20:23

  A <- 0.0002
  B <- 0.1
  C <- 0.3
  D <- 0.03
  E <- 10
  F_hp <- 20
  G <- 0.00005
  H <- 1.08

  tab <- mortality_law_table(
    "HeligmanPollard",
    min(ages),
    max(ages),
    A = A,
    B = B,
    C = C,
    D = D,
    E = E,
    F_hp = F_hp,
    G = G,
    H = H,
    close = FALSE
  )

  xx <- as.numeric(ages)

  odds <- (A^((xx + B)^C)) +
    (D * exp(-E * (log(xx) - log(F_hp))^2)) +
    (G * (H^xx))

  expected_q <- odds / (1 + odds)

  expect_true(all(is.na(tab$mu_x)))
  expect_equal(tab$qx, expected_q, tolerance = 1e-12)
})

test_that("Heligman-Pollard transitional F alias agrees with F_hp", {
  args <- list(
    law = "HeligmanPollard",
    x_min = 20,
    x_max = 23,
    A = 0.0002,
    B = 0.1,
    C = 0.3,
    D = 0.03,
    E = 10,
    G = 0.00005,
    H = 1.08,
    close = FALSE
  )

  current <- do.call(
    mortality_law_table,
    c(args, list(F_hp = 20))
  )

  legacy <- do.call(
    mortality_law_table,
    c(args, list(F = 20))
  )

  expect_equal(legacy$qx, current$qx, tolerance = 1e-12)
})

test_that("mortality_law_table closure and life-table accounting are coherent", {
  tab <- mortality_law_table(
    "Exponential",
    60,
    64,
    lambda = 0.02,
    l0 = 100000,
    a_x = 0.5,
    close = TRUE
  )

  expect_equal(tail(tab$qx, 1), 1)
  expect_equal(tail(tab$px, 1), 0)

  expect_equal(
    tab$dx,
    tab$lx * tab$qx,
    tolerance = 1e-10
  )

  lx_next <- c(tab$lx[-1], tab$lx[length(tab$lx)] * tab$px[length(tab$px)])

  expect_equal(
    tab$Lx,
    lx_next + 0.5 * tab$dx,
    tolerance = 1e-10
  )

  expect_equal(
    tab$Tx,
    rev(cumsum(rev(tab$Lx))),
    tolerance = 1e-10
  )

  expect_equal(
    tab$ex,
    tab$Tx / tab$lx,
    tolerance = 1e-10
  )
})

test_that("params can be supplied as a named list and dots override them", {
  from_params <- mortality_law_table(
    "Gompertz",
    40,
    43,
    params = list(B = 1e-5, c = 1.08),
    close = FALSE
  )

  overridden <- mortality_law_table(
    "Gompertz",
    40,
    43,
    c = 1.10,
    params = list(B = 1e-5, c = 1.08),
    close = FALSE
  )

  expect_equal(
    from_params$mu_x,
    1e-5 * 1.08^(40:43),
    tolerance = 1e-12
  )

  expect_equal(
    overridden$mu_x,
    1e-5 * 1.10^(40:43),
    tolerance = 1e-12
  )
})

test_that("mortality_law_table validates essential law parameters", {
  expect_error(
    mortality_law_table(
      "Gompertz",
      40,
      50,
      B = 1e-5
    ),
    "Missing parameter"
  )

  expect_error(
    mortality_law_table(
      "Exponential",
      40,
      50,
      lambda = -0.01
    ),
    "lambda > 0"
  )

  expect_error(
    mortality_law_table(
      "DeMoivre",
      60,
      70,
      omega = 65
    ),
    "omega >= x_max"
  )

  expect_error(
    mortality_law_table(
      "Beta",
      0,
      100,
      alpha = 2,
      beta = 5,
      omega = 100
    ),
    "omega >= x_max \\+ 1"
  )
})
