test_that("single-life values retain original rates and representable products", {
  lt <- data.frame(x = 60:62, lx = c(1, 1e-300, 0))
  a <- annuity_x(lt, 60, i = -1000, i_type = "force", n = 1)
  expect_equal(a / exp(log(1e-300) + 1000), 1)
  expect_equal(annuity_x(lt, 60, i = 1000, i_type = "force",
                         n = 1, timing = "due"), 1)
  lt <- data.frame(x = 60:62, lx = c(100, 50, 0))
  b <- insurance_x(lt, 60, i = 800, i_type = "force", n = 1,
                   type = "term", benefit = 1e300)
  expect_equal(b / exp(log(1e300) + log(.5) - 800), 1)
  b <- insurance_x(lt, 60, i = -800, i_type = "force", n = 1,
                   type = "term", benefit = -1e-300)
  expect_equal(b / -exp(log(1e-300) + log(.5) + 800), 1)
  expect_equal(insurance_x(lt, 60, i = -800, i_type = "force",
                           n = 1, type = "term", benefit = 0), 0)
})

test_that("deferred values combine survival and discount before rounding", {
  lt <- data.frame(x = 60:62, lx = c(1e300, 1e-100, 0))
  a <- annuity_x(lt, 60, i = -800, i_type = "force", h = 1,
                 n = 1, timing = "due", tidy = TRUE)
  expect_equal(a$apv / exp(log(1e-100) - log(1e300) + 800), 1)
  b <- insurance_x(lt, 60, i = -800, i_type = "force", h = 1,
                   n = 1, type = "term", benefit = 1e-200, tidy = TRUE)
  ref <- exp(log(1e-200) + log(1e-100) - log(1e300) + 1600)
  expect_equal(b$apv / ref, 1)
})

test_that("extinct tails contribute zero instead of undefined conditional deaths", {
  lt <- data.frame(x = 60:63, lx = c(100, 0, 0, 0))
  expect_equal(insurance_x(lt, 60, i = .05), 1 / 1.05)
  expect_equal(insurance_x(lt, 60, i = .05, h = 1), 0)
  for (frac in c("UDD", "CF", "CML", "Balducci")) {
    a <- annuity_x(lt, 60, i = .05, n = 3, k = 4, frac = frac)
    p <- if (frac == "UDD") c(.75, .5, .25) else rep(0, 3)
    expect_equal(a, sum(p * 1.05^-(1:3 / 4)) / 4)
    expect_equal(annuity_x(lt, 60, i = .05, h = 1, n = 1,
                           k = 4, frac = frac), 0)
  }
})

test_that("annual insurance uses survivor decrements without ratio cancellation", {
  lt <- data.frame(x = 60:62, lx = c(1e16, 1e16 - 2, 0))
  observed <- insurance_x(lt, 60, i = .05, n = 1, type = "term")
  expect_equal(observed / (2e-16 / 1.05), 1, tolerance = 1e-13)
})

test_that("exact single-life values match independent payment enumerations", {
  lt <- data.frame(x = 60:65, lx = c(1000, 800, 500, 200, 50, 0))
  for (h in 0:1) for (i in c(-.03, 0, .05)) {
    n <- 2
    death <- (lt$lx[h + 1:2] - lt$lx[h + 2:3]) / 1000
    survival <- lt$lx[h + 1:3] / 1000
    expected_a <- sum(survival[1:2] * (1 + i)^-(h + 0:1))
    expected_b <- sum(death * (1 + i)^-(h + 1:2))
    endowment <- survival[3] * (1 + i)^-(h + n)
    a <- annuity_x(lt, 60, i, n = n, h = h, timing = "due")
    b <- insurance_x(lt, 60, i, n = n, h = h, type = "term")
    expect_equal(a, expected_a)
    expect_equal(b, expected_b)
    expect_equal(insurance_x(lt, 60, i, n = n, h = h,
                             type = "endowment"), expected_b + endowment)
    expect_equal(b + i / (1 + i) * a,
                 survival[1] * (1 + i)^-h - endowment)
  }
  for (frac in c("UDD", "CF", "CML", "Balducci")) {
    times <- (1:6) / 4
    j <- floor(times)
    s <- times - j
    p <- lt$lx[j + 2] / lt$lx[j + 1]
    within <- switch(frac, UDD = 1 - s * (1 - p), CF = p^s,
                     CML = p^s, Balducci = p / (p + s * (1 - p)))
    survival <- lt$lx[j + 1] / 1000 * within
    ref <- sum(survival * 1.04^-times) / 4
    expect_equal(annuity_x(lt, 60, .04, n = 1.5, k = 4, frac = frac), ref)
  }
})

test_that("Woolhouse orders match the independent constant-force expansion", {
  mu <- .03
  delta <- .04
  lt <- data.frame(x = 60:80, lx = 1000 * exp(-mu * (0:20)))
  n <- 3
  k <- 12
  e <- exp(-(mu + delta) * n)
  annual <- -expm1(-(mu + delta) * n) / -expm1(-(mu + delta))
  first <- annual - (k - 1) / (2 * k) * (1 - e)
  second <- first - (k^2 - 1) / (12 * k^2) * (mu + delta) * (1 - e)
  for (order in c("first", "second")) {
    ref <- if (order == "first") first else second
    expect_equal(annuity_x(lt, 60, delta, i_type = "force", n = n, k = k,
                           timing = "due", woolhouse = order), ref)
    expect_equal(annuity_x(lt, 60, delta, i_type = "force", n = n, k = k,
                           timing = "immediate", woolhouse = order), ref - (1 - e) / k)
  }
})

test_that("insurance respects explicitly supplied contract interest overrides", {
  lt <- data.frame(x = 60:62, lx = c(100, 50, 0))
  contract <- life_contract(lt, x = 60, i = .2, i_type = "force")
  expect_equal(insurance_x(contract, i = .1, i_type = "effective", n = 1,
                           type = "term"), .5 / 1.1)
})

test_that("single-life inputs fail clearly before integer coercion or valuation", {
  lt <- data.frame(x = 60:62, lx = c(100, 80, 0))
  for (fun in list(annuity_x, insurance_x)) {
    expect_error(fun(lt, x = 1e20, i = .05), "x")
    expect_error(fun(lt, x = 59, i = .05), "age|x")
    expect_error(fun(lt, x = 62, i = .05), "zero|positive|lx")
    expect_error(fun(transform(lt, lx = c(100, 110, 0)), 60, .05), "nonincreasing|non-increasing")
    expect_error(fun(lt[FALSE, ], 60, .05), "empty|age|row")
    expect_no_warning(out <- fun(lt, 60, .05, m = 3e9))
    expect_true(is.finite(out))
  }
  expect_error(annuity_x(lt, 60, .05, n = 1, k = 3e9), "payments|payment")
})

test_that("terminal ages and truly zero payments do not overflow intermediate arithmetic", {
  lt <- data.frame(x = c(2147483646, 2147483647), lx = c(100, 50))
  expect_equal(insurance_x(lt, 2147483646, 0), 1)
  expect_equal(annuity_x(lt, 2147483646, 0, timing = "due"), 1.5)
  lt <- data.frame(x = 60:63, lx = c(100, 0, 0, 0))
  expect_equal(annuity_x(lt, 60, -1e308, i_type = "force", timing = "due"), 1)
  expect_equal(insurance_x(lt, 60, -1e308, i_type = "force", benefit = 0), 0)
  expect_equal(insurance_x(lt, 60, -1e308, i_type = "force", h = 2), 0)
})

test_that("Woolhouse reports failures of its smooth positive-survival approximation", {
  lt <- data.frame(x = 60:65, lx = 1000 * exp(-.03 * (0:5)))
  expect_warning(value <- annuity_x(lt, 60, 10, i_type = "force", n = 1,
                                     k = 12, timing = "due", woolhouse = "second"),
                 "negative")
  expect_lt(value, 0) # Report the approximation; do not silently clamp its value.
  lt <- data.frame(x = 60:63, lx = c(100, 0, 0, 0))
  expect_error(annuity_x(lt, 60, .05, n = 1, k = 12, woolhouse = "second"),
               "positive.*survival")
})
