test_that("forward_rate satisfies the actuarial spot-forward identity", {
  t <- c(1, 2, 3, 5)
  i <- c(0.035, 0.04, 0.047, 0.055)
  t_start <- 2
  t_end <- 5

  expected <- (
    (1 + i[4])^t_end /
      (1 + i[2])^t_start
  )^(1 / (t_end - t_start)) - 1

  out <- forward_rate(
    t = t,
    i = i,
    t_start = t_start,
    t_end = t_end,
    method = "exact"
  )

  expect_equal(out$f, expected, tolerance = 1e-12)
  expect_equal(out$i_start, i[2], tolerance = 1e-12)
  expect_equal(out$i_end, i[4], tolerance = 1e-12)

  lhs <- (1 + out$i_start)^t_start *
    (1 + out$f)^(t_end - t_start)
  rhs <- (1 + out$i_end)^t_end

  expect_equal(lhs, rhs, tolerance = 1e-12)
})


test_that("forward_rate standardizes nominal spot rates before calculation", {
  t <- c(1, 2, 3)
  j <- c(0.05, 0.055, 0.06)
  m <- 2
  t_start <- 1
  t_end <- 3

  i_eff <- (1 + j / m)^m - 1
  expected <- (
    (1 + i_eff[3])^3 /
      (1 + i_eff[1])^1
  )^(1 / 2) - 1

  out <- forward_rate(
    t = t,
    i = j,
    t_start = t_start,
    t_end = t_end,
    i_type = "nominal_interest",
    m = m
  )

  expect_equal(out$i_start, i_eff[1], tolerance = 1e-12)
  expect_equal(out$i_end, i_eff[3], tolerance = 1e-12)
  expect_equal(out$f, expected, tolerance = 1e-12)
})


test_that("forward_rate linear method interpolates spot rates before forming the forward", {
  t <- c(1, 3, 5)
  i <- c(0.03, 0.05, 0.07)
  t_start <- 2
  t_end <- 4

  i_start <- 0.04
  i_end <- 0.06
  expected <- (
    (1 + i_end)^t_end /
      (1 + i_start)^t_start
  )^(1 / (t_end - t_start)) - 1

  out <- forward_rate(
    t = t,
    i = i,
    t_start = t_start,
    t_end = t_end,
    method = "linear"
  )

  expect_equal(out$i_start, i_start, tolerance = 1e-12)
  expect_equal(out$i_end, i_end, tolerance = 1e-12)
  expect_equal(out$f, expected, tolerance = 1e-12)
})


test_that("forward_rate exact method requires observed maturity nodes", {
  expect_error(
    forward_rate(
      t = c(1, 3, 5),
      i = c(0.03, 0.05, 0.07),
      t_start = 2,
      t_end = 5,
      method = "exact"
    ),
    "does not match any curve node"
  )
})


test_that("forward_rate never extrapolates beyond the observed curve", {
  expect_error(
    forward_rate(
      t = c(1, 3, 5),
      i = c(0.03, 0.05, 0.07),
      t_start = 2,
      t_end = 6,
      method = "linear"
    ),
    "outside the observed maturity range"
  )
})


test_that("forward_rate handles multiple curves in a tidy workflow", {
  curves <- tibble::tibble(
    curve = c("base", "stress"),
    t = list(c(1, 2, 4), c(1, 3, 5)),
    i = list(c(0.03, 0.04, 0.05), c(0.05, 0.055, 0.06)),
    t_start = c(1, 1),
    t_end = c(4, 5)
  )

  out <- forward_rate(curves)

  expected_1 <- ((1.05^4) / (1.03^1))^(1 / 3) - 1
  expected_2 <- ((1.06^5) / (1.05^1))^(1 / 4) - 1

  expect_equal(nrow(out), 2L)
  expect_equal(out$curve, c("base", "stress"))
  expect_equal(out$f, c(expected_1, expected_2), tolerance = 1e-12)
})


test_that("forward_rate supports custom columns and computed-only output", {
  curves <- tibble::tibble(
    maturities = list(c(1, 2, 3)),
    spots = list(c(0.04, 0.045, 0.05)),
    start = 1,
    end = 3
  )

  out <- forward_rate(
    curves,
    col_t = "maturities",
    col_i = "spots",
    col_t_start = "start",
    col_t_end = "end",
    .out = "forward",
    .keep = "none"
  )

  expected <- ((1.05^3) / 1.04)^(1 / 2) - 1

  expect_named(out, c("forward", "i_start", "i_end"))
  expect_equal(out$forward, expected, tolerance = 1e-12)
})


test_that("forward_rate plot output is a ggplot object", {
  out <- forward_rate(
    t = c(1, 2, 3),
    i = c(0.04, 0.045, 0.05),
    t_start = 1,
    t_end = 3,
    plot = TRUE
  )

  expect_true("forward_rate_plot" %in% names(out))
  expect_s3_class(out$forward_rate_plot[[1]], "ggplot")
})


test_that("forward_rate propagates missing values by row", {
  curves <- tibble::tibble(
    curve = c("complete", "missing"),
    t = list(c(1, 2, 3), c(1, 2, 3)),
    i = list(c(0.04, 0.045, 0.05), c(0.04, NA_real_, 0.05)),
    t_start = c(1, 1),
    t_end = c(3, 3)
  )

  out <- forward_rate(curves, .na = "propagate")

  expect_false(is.na(out$f[1]))
  expect_true(is.na(out$f[2]))
})


test_that("forward_rate validates the ordering of the forward interval", {
  expect_error(
    forward_rate(
      t = c(1, 2, 3),
      i = c(0.04, 0.045, 0.05),
      t_start = 2,
      t_end = 2
    ),
    "strictly greater"
  )
})
