test_that("nearby interest and growth match explicit payments", {
  r <- .03
  for (s in r + c(-1e-10, -5e-13, 0, 5e-13, 1e-10)) {
    for (timing in c("immediate", "due")) {
      j <- 0:47
      times <- (j + as.integer(timing == "immediate")) / 12
      pv <- sum(100 * exp(s * j / 12 - r * (2 + times)))
      fv <- sum(100 * exp(s * j / 12 + r * (4 - times)))
      out <- annuity_geom(4, k = 12, i = r, i_type = "force", g = s,
        g_type = "force", h = 2, P1 = 100, timing = timing, tidy = TRUE)
      expect_equal(out$present_value_factor, pv, tolerance = 1e-12)
      expect_equal(out$accumulated_value_factor, fv, tolerance = 1e-12)
    }
  }
})

test_that("a small growth difference is retained over a long term", {
  # Decimal 80-digit integral references; discrete-relative error <=exp(5e-13)-1.
  for (s in c(-5e-13, 5e-13)) {
    reference <- if (s < 0) 786938680574.733152792400930017639093 else
      1297442541400.25629369730157562832714
    for (valuation in c("present", "accumulated")) {
      expect_equal(annuity_geom(1e12, i = 0, g = s, g_type = "force",
        valuation = valuation) / reference, 1, tolerance = 1e-12)
    }
  }
  # Decimal on exact inputs: delta_i=1, delta_g=1+2^-52, k=10, N=1e17.
  # Subtract annual log rates before division; subtracting rounded period
  # rates changes the exponent N*(s-r) materially in this example.
  expect_equal(annuity_geom(1e16, k = 10, i = 1, i_type = "force",
    g = 1 + 2^-52, g_type = "force") /
    334618217753837351.97066129644432654, 1, tolerance = 1e-12)
})

test_that("original forces and monetary scale survive overflowing intermediates", {
  force <- 1100 * log(2)
  expect_equal(annuity_geom(2, i = force, i_type = "force", g = force,
    g_type = "force", timing = "due"), 2)
  expect_equal(annuity_geom(2, i = -force, i_type = "force", g = -force,
    g_type = "force", P1 = 2^-1000) / 2^101, 1, tolerance = 1e-12)
  expect_equal(annuity_geom(2, i = -force, i_type = "force", g = -force,
    g_type = "force", P1 = 2^1000, valuation = "accumulated") /
    2^-99, 1, tolerance = 1e-12)
  expect_equal(annuity_geom(2, i = 800, i_type = "force", g = -800,
    g_type = "force", P1 = 0, tidy = TRUE)$accumulated_value_factor, 0)
  expect_equal(annuity_geom(1, i = 800, i_type = "force", g = -800,
    g_type = "force", timing = "due"), 1)
})

test_that("perpetuity convergence uses strict log-rate order", {
  expect_equal(annuity_geom(i = 2^-42, g = 0, perpetuity = TRUE) / 2^42,
    1, tolerance = 1e-12)
  expect_equal(annuity_geom(i = -.01, g = -.02, perpetuity = TRUE),
    100, tolerance = 1e-12)
  expect_equal(annuity_geom(i = 800, i_type = "force", g = 799,
    g_type = "force", perpetuity = TRUE, timing = "due"),
    1 / (1 - exp(-1)), tolerance = 1e-12)
  expect_error(annuity_geom(i = .03, g = .03, perpetuity = TRUE), "strictly greater")
  expect_error(annuity_geom(i = .03, g = .04, perpetuity = TRUE), "strictly greater")
})

test_that("geometric scenarios preserve missing values and editable table workflows", {
  tab <- tibble::tibble(term = c(2, 3), first = c(10, -20))
  out <- tab |> dplyr::mutate(value = annuity_geom(term, i = .04, g = .02, P1 = first))
  expected <- c(sum(10 * 1.02^(0:1) / 1.04^(1:2)),
    sum(-20 * 1.02^(0:2) / 1.04^(1:3)))
  expect_equal(out$value, expected, tolerance = 1e-12)
  expect_true(is.na(annuity_geom(NA_real_, i = .04)))
  expect_true(is.na(annuity_geom(2, i = .04, g_type = NA_character_)))
  expect_true(is.na(annuity_geom(2, i = .04, P1 = NA_real_)))
})

test_that("nominal growth and interest use separate conversion frequencies", {
  delta_i <- 12 * log1p(.06 / 12)
  delta_g <- -2 * log1p(-.024 / 2)
  times <- (1:12) / 4
  out <- annuity_geom(3, k = 4, i = .06, i_type = "nominal_interest", m = 12,
    g = .024, g_type = "nominal_discount", g_m = 2, h = .5, P1 = -25,
    valuation = "accumulated", tidy = TRUE)
  expect_equal(out$present_value_factor,
    sum(-25 * exp(delta_g * (0:11) / 4 - delta_i * (.5 + times))), tolerance = 1e-12)
  expect_equal(out$annuity_factor,
    sum(-25 * exp(delta_g * (0:11) / 4 + delta_i * (3 - times))), tolerance = 1e-12)
  expect_error(annuity_geom(2, i = .04, m = 0), "`m`")
  expect_error(annuity_geom(2, i = .04, g_m = 0), "`g_m`")
})

test_that("matched rate families retain close spreads before logarithms", {
  i <- .06
  g <- i - 1e-12
  expect_equal(annuity_geom(i = i, g = g, perpetuity = TRUE) /
    (1 / (i - g)), 1, tolerance = 1e-12)
  expect_equal(annuity_geom(k = 12, i = i, i_type = "nominal_interest", m = 12,
    g = g, g_type = "nominal_interest", g_m = 12, perpetuity = TRUE) /
    (12 / (i - g)), 1, tolerance = 1e-12)
  expected <- (12 - i) * (12 - g) / (12 * (i - g))
  expect_equal(annuity_geom(k = 12, i = i, i_type = "nominal_discount", m = 12,
    g = g, g_type = "nominal_discount", g_m = 12, perpetuity = TRUE) /
    expected, 1, tolerance = 1e-12)
  # Decimal 80 digits: exact effective i=1/2, g=1/2+2^-52, N=1e16.
  expect_equal(annuity_geom(1e16, i = .5, g = .5 + 2^-52) /
    15286352891704535.0090512702373029173, 1, tolerance = 1e-12)
  first <- 2^-1000
  # A widely separated ratio overflows; fall back to separate finite logs.
  expect_equal(annuity_geom(2, i = -.5, g = 1e308, P1 = first) /
    (2 * first + 4 * (first * 1e308)), 1, tolerance = 1e-12)
  # A nominal ratio denominator can itself overflow; it must not erase a
  # large negative log spread. Only the first due payment then contributes.
  expect_equal(annuity_geom(2, i = 1e308, i_type = "nominal_interest", m = 1e308,
    g = 1e308 - 2e292, g_type = "nominal_interest", g_m = 1e308,
    timing = "due"), 1)
  expect_equal(annuity_geom(2, i = -1e308, i_type = "nominal_discount", m = 1e308,
    g = -1e308 - 2e292, g_type = "nominal_discount", g_m = 1e308,
    timing = "due"), 1)
})
