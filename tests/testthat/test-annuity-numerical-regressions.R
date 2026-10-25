test_that("small rates retain their effect over long terms", {
  # Independent continuous reference: integrate exp(-u/2) on [0,1],
  # then rescale time by 1e12. Decimal references are recorded in the audit.
  n <- 1e12
  delta <- 5e-13
  pv_cont <- 786938680574.733152792400930017639093
  fv_cont <- 1297442541400.25629369730157562832714
  for (sgn in c(-1, 1)) {
    pv <- if (sgn > 0) pv_cont else fv_cont
    fv <- if (sgn > 0) fv_cont else pv_cont
    expect_equal(a_angle(n, i = sgn * delta, i_type = "force", timing = "continuous") /
      pv, 1, tolerance = 1e-12)
    expect_equal(s_angle(n, i = sgn * delta, i_type = "force", timing = "continuous") /
      fv, 1, tolerance = 1e-12)
    # A Riemann sum with mesh 1e-12 differs relatively from this integral
    # by at most exp(abs(delta))-1, below the tolerance used here.
    expect_equal(a_angle(n, i = sgn * delta, i_type = "force") / pv,
      1, tolerance = 1e-12)
    expect_equal(s_angle(n, i = sgn * delta, i_type = "force") / fv,
      1, tolerance = 1e-12)
  }
})

test_that("vectorized tabular workflows preserve NA and level-payment units", {
  cases <- tibble::tibble(term = c(2, 3), frequency = c(2, 12),
    force = c(-.02, .04), payment = c(100, 250))
  values <- cases |> dplyr::mutate(
    factor = a_angle(term, k = frequency, i = force, i_type = "force"),
    value = payment * factor)
  expected <- vapply(seq_len(nrow(cases)), function(j) {
    times <- seq_len(cases$term[j] * cases$frequency[j]) / cases$frequency[j]
    sum(cases$payment[j] * exp(-cases$force[j] * times))
  }, numeric(1))
  expect_s3_class(values, "tbl_df")
  expect_equal(values$value, expected, tolerance = 1e-12)
  for (fun in list(a_angle, s_angle)) {
    expect_equal(c(2, 3) |> fun(i = 0), c(2, 3))
    expect_true(is.na(fun(2, i = .03, i_type = NA_character_)))
    out <- fun(2, i = NA_real_, payment = 0, tidy = TRUE)
    value <- if ("present_value" %in% names(out)) out$present_value else out$future_value
    expect_true(is.na(value))
  }
  expect_error(a_angle(i = 0, perpetuity = TRUE), "strictly positive")
  expect_error(a_angle(i = -.01, perpetuity = TRUE, timing = "continuous"),
    "strictly positive")
})

test_that("ordinary annuities match explicit payments near zero", {
  for (delta in c(-1e-10, 0, 1e-10, 0.03)) {
    for (timing in c("immediate", "due")) {
      k <- 12
      times <- (seq_len(36) - as.integer(timing == "due")) / k
      expect_equal(a_angle(3, k = k, i = delta, i_type = "force", h = 2,
        timing = timing), sum(exp(-delta * (2 + times))), tolerance = 1e-12)
      expect_equal(s_angle(3, k = k, i = delta, i_type = "force", h = 2,
        timing = timing), sum(exp(delta * (3 - times))), tolerance = 1e-12)
    }
  }
})

test_that("finite original forces survive saturated annual equivalents", {
  for (delta in c(-800, 800)) {
    # Exactly one due payment at issue, and one immediate payment at horizon.
    expect_equal(a_angle(1, i = delta, i_type = "force", timing = "due"), 1)
    expect_equal(s_angle(1, i = delta, i_type = "force"), 1)
    expect_equal(a_angle(0.5, k = 1600, i = delta, i_type = "force", h = 0,
      timing = if (delta > 0) "immediate" else "due") /
      sum(exp(-delta * ((seq_len(800) - as.integer(delta < 0)) / 1600))),
      1, tolerance = 1e-12)
  }
  # Known zero monetary value with a genuinely overflowing annuity factor.
  expect_equal(s_angle(2, i = 800, i_type = "force", payment = 0,
    tidy = TRUE)$future_value, 0)
})

test_that("positive perpetuities have no arbitrary lower rate cutoff", {
  delta <- 5e-13
  expect_equal(a_angle(i = delta, i_type = "force", perpetuity = TRUE,
    timing = "continuous") / 2e12, 1, tolerance = 1e-12)
  # sum exp(-delta*j), j>=1, is 1/expm1(delta).
  expect_equal(a_angle(i = delta, i_type = "force", perpetuity = TRUE) /
    (1 / expm1(delta)), 1, tolerance = 1e-12)
  expect_equal(a_angle(i = 800, i_type = "force", perpetuity = TRUE,
    timing = "due"), 1)
})

test_that("continuous integrals retain tiny terms and bounded large spans", {
  for (fun in list(a_angle, s_angle)) {
    for (delta in c(-1e-300, 1e-300)) {
      expect_equal(fun(1e-20, i = delta, i_type = "force", timing = "continuous") /
        1e-20, 1, tolerance = 1e-12)
    }
  }
  # Integrals of exp(-1e308*t) have total 1e-308; omitted tails are negligible.
  expect_equal(a_angle(1e308, i = 1e308, i_type = "force", timing = "continuous") /
    1e-308, 1, tolerance = 1e-12)
  expect_equal(s_angle(1e308, i = -1e308, i_type = "force", timing = "continuous") /
    1e-308, 1, tolerance = 1e-12)
})

test_that("monetary values can remain finite when unit factors lose range", {
  force <- 1100 * log(2)
  expect_equal(a_angle(1, i = -force, i_type = "force", payment = 2^-1000,
    tidy = TRUE)$present_value / 2^100, 1, tolerance = 1e-12)
  expect_equal(a_angle(1, i = force, i_type = "force", payment = 2^1000,
    tidy = TRUE)$present_value / 2^-100, 1, tolerance = 1e-12)
  expect_equal(s_angle(2, i = force, i_type = "force", payment = 2^-1000,
    tidy = TRUE)$future_value / 2^100, 1, tolerance = 1e-12)
  expect_equal(s_angle(1, i = -force, i_type = "force", timing = "due",
    payment = 2^1000, tidy = TRUE)$future_value / 2^-100,
    1, tolerance = 1e-12)
})
