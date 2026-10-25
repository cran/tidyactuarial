test_that("forward_curve recovers each adjacent effective forward rate", {
  t <- c(1, 2, 3)
  i <- c(0.04, 0.045, 0.05)
  v <- (1 + i)^(-t)
  expected <- c(
    0.04,
    (v[1] / v[2])^(1 / (t[2] - t[1])) - 1,
    (v[2] / v[3])^(1 / (t[3] - t[2])) - 1
  )
  out <- forward_curve(t = t, i = i)
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("curve", "t_start", "t_end", "duration",
                      "v_start", "v_end", "f"))
  expect_equal(out$curve, rep("1", 3))
  expect_equal(out$t_start, c(0, 1, 2))
  expect_equal(out$t_end, t)
  expect_equal(out$f, expected, tolerance = 1e-12)
})

test_that("forward_curve reprices the original spot discount factors", {
  t <- c(0.5, 1, 2.5, 4)
  i <- c(0.025, 0.03, 0.037, 0.045)
  out <- forward_curve(t = t, i = i)
  implied <- cumprod((1 + out$f)^(-out$duration))
  expect_equal(implied, (1 + i)^(-t), tolerance = 1e-12)
  expect_equal(out$duration, diff(c(0, t)))
})

test_that("forward_curve handles constant spot curves with irregular maturities", {
  t <- c(0.25, 1.5, 2, 7)
  out <- forward_curve(t = t, i = rep(0.052, length(t)))
  expect_equal(out$f, rep(0.052, length(t)), tolerance = 1e-12)
})

test_that("forward_curve accepts yield_curve outputs with list columns", {
  curve <- yield_curve(t = 1:3, i = c(0.04, 0.045, 0.05))
  actual <- forward_curve(curve)
  reference <- forward_curve(t = 1:3, i = c(0.04, 0.045, 0.05))
  expect_equal(actual, reference, tolerance = 1e-12)
})

test_that("forward_curve accepts long curves including bootstrap output", {
  t <- 1:3
  i <- c(0.04, 0.045, 0.05)
  v <- (1 + i)^(-t)
  bonds <- tibble::tibble(
    instrument = c("Z1", "B2", "B3"),
    price = c(100 * v[1], 5 * v[1] + 105 * v[2],
              6 * v[1] + 6 * v[2] + 106 * v[3]),
    t = list(1, c(1, 2), c(1, 2, 3)),
    cf = list(100, c(5, 105), c(6, 6, 106))
  )
  actual <- bootstrap_spot_curve(bonds) |>
    forward_curve()
  reference <- forward_curve(t = t, i = i)
  expect_equal(actual, reference, tolerance = 1e-11)
})

test_that("forward_curve sorts unsorted numeric nodes and custom column names", {
  x <- tibble::tibble(maturity = c(3, 1, 2),
                      spot = c(0.05, 0.04, 0.045))
  original <- x
  out <- forward_curve(x, col_t = "maturity", col_i = "spot")
  expect_equal(out$t_end, c(1, 2, 3))
  expect_equal(out$f,
               forward_curve(t = 1:3, i = c(0.04, 0.045, 0.05))$f,
               tolerance = 1e-12)
  expect_identical(x, original)
})

test_that("forward_curve separates curves by a grouping column", {
  x <- tibble::tibble(
    scenario = c("base", "stress", "base", "stress"),
    maturity = c(2, 1, 1, 2),
    spot = c(0.05, 0.025, 0.04, 0.03)
  )
  out <- forward_curve(x, col_t = "maturity", col_i = "spot",
                       col_curve = "scenario")
  expect_equal(out$curve, c("base", "base", "stress", "stress"))
  expect_equal(out$t_end, c(1, 2, 1, 2))
  expect_equal(out$f[1:2],
               forward_curve(t = 1:2, i = c(0.04, 0.05))$f,
               tolerance = 1e-12)
})

test_that("forward_curve preserves identifiers for multiple list-column curves", {
  x <- tibble::tibble(
    label = c("base", "stress"),
    t = list(c(1, 2, 3), c(1, 3, 5)),
    i = list(c(0.04, 0.045, 0.05), c(0.05, 0.055, 0.06))
  )
  out <- forward_curve(x, col_curve = "label")
  expect_equal(out$curve, c(rep("base", 3), rep("stress", 3)))
  expect_equal(out$t_start[4:6], c(0, 1, 3))
  expect_equal(out$t_end[4:6], c(1, 3, 5))
})

test_that("forward_curve handles explicit zero maturity and suppresses origin", {
  with_zero <- forward_curve(t = c(0, 1, 2), i = c(0, 0.04, 0.05))
  without_zero <- forward_curve(t = c(1, 2), i = c(0.04, 0.05))
  expect_equal(with_zero, without_zero, tolerance = 1e-12)
  no_origin <- forward_curve(t = 1:3, i = c(0.04, 0.045, 0.05),
                             include_origin = FALSE)
  expect_equal(no_origin$t_start, c(1, 2))
  expect_equal(no_origin$t_end, c(2, 3))
  expect_equal(no_origin$f,
               forward_curve(t = 1:3, i = c(0.04, 0.045, 0.05))$f[2:3],
               tolerance = 1e-12)
})

test_that("forward_curve supports nominal interest and force of interest", {
  annual <- forward_curve(t = c(1, 2), i = c(0.06, 0.08))
  nominal <- forward_curve(t = c(1, 2),
                           i = 2 * (sqrt(1 + c(0.06, 0.08)) - 1),
                           i_type = "nominal_interest", m = 2)
  force <- forward_curve(t = c(1, 2), i = log(1 + c(0.06, 0.08)),
                         i_type = "force")
  expect_equal(nominal$f, annual$f, tolerance = 1e-12)
  expect_equal(force$f, annual$f, tolerance = 1e-12)
})

test_that("forward_curve supports negative spot and forward rates", {
  x <- forward_curve(t = 1:3, i = c(-0.02, -0.03, -0.04))
  expect_true(all(x$f < 0))
  expect_true(all(x$f > -1))
  expect_equal(cumprod((1 + x$f)^(-x$duration)),
               (1 + c(-0.02, -0.03, -0.04))^(-(1:3)),
               tolerance = 1e-12)
})

test_that("forward_curve is consistent with existing forward_rate", {
  t <- c(1, 2, 3)
  i <- c(0.04, 0.045, 0.05)
  curve <- forward_curve(t = t, i = i)
  selected <- forward_rate(t = t, i = i, t_start = 1, t_end = 2)
  expect_equal(curve$f[2], selected$f, tolerance = 1e-12)
})

test_that("forward_curve rejects invalid maturities and spot rates", {
  expect_error(forward_curve(t = c(1, 1), i = c(0.04, 0.05)), "duplicate")
  expect_error(forward_curve(t = c(1, NA), i = c(0.04, 0.05)), "finite")
  expect_error(forward_curve(t = c(1, Inf), i = c(0.04, 0.05)), "finite")
  expect_error(forward_curve(t = c(-1, 2), i = c(0.04, 0.05)), "nonnegative")
  expect_error(forward_curve(t = 1:2, i = c(0.04, NA)), "finite")
  expect_error(forward_curve(t = 1:2, i = c(0.04, -1)), "greater than -1")
  expect_error(forward_curve(t = 1:2, i = 0.04), "equal-length")
})

test_that("forward_curve validates curve inputs and missing identifiers", {
  expect_error(forward_curve(t = 1, i = 0.04, include_origin = FALSE),
               "fewer than two nodes")
  expect_error(forward_curve(t = 0, i = 0.04), "fewer than two nodes")
  expect_error(forward_curve(t = 1:2, i = c(0.04, 0.05),
                             include_origin = NA), "TRUE or FALSE")
  expect_error(forward_curve(t = 1:2, i = c(0.04, 0.05),
                             col_curve = "scenario"), "requires")
  expect_error(forward_curve(tibble::tibble(t = 1:2, i = c(0.04, 0.05)),
                             col_curve = "scenario"), "Missing required")
  expect_error(forward_curve(tibble::tibble(
    scenario = c("A", NA_character_), t = 1:2,
    i = c(0.04, 0.05)), col_curve = "scenario"), "nonmissing")
})

test_that("forward_curve refuses ambiguous list-curve identifiers", {
  x <- tibble::tibble(
    scenario = c("same", "same"),
    t = list(1:2, 1:2),
    i = list(c(0.04, 0.05), c(0.03, 0.04))
  )
  expect_error(forward_curve(x, col_curve = "scenario"), "unique")
  expect_error(forward_curve(x, i_type = c("effective", "force")), "scalars")
})

test_that("forward_curve rejects impossible discount factor calculations", {
  expect_error(forward_curve(t = c(1, 10000), i = c(0.04, 0.30)),
               "nonpositive or nonfinite")
  expect_error(forward_curve(t = c(1, 2), i = c(0.04, 0.05), m = 0),
               "positive integer")
  expect_error(forward_curve(t = c(1, 2), i = c(0.04, 0.05),
                             i_type = "unknown"), "i_type")
})
