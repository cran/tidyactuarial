test_that("yield_curve builds spot discount factors from effective spot rates", {
  t <- c(1, 2, 3, 5)
  i <- c(0.04, 0.045, 0.05, 0.055)

  out <- yield_curve(t = t, i = i)
  expected <- (1 + i)^(-t)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
  expect_equal(out$t[[1]], t)
  expect_equal(out$i[[1]], i)
  expect_equal(out$i_effective[[1]], i, tolerance = 1e-12)
  expect_equal(out$v[[1]], expected, tolerance = 1e-12)
})


test_that("yield_curve standardizes nominal spot rates before discounting", {
  t <- c(1, 2, 4)
  j <- c(0.06, 0.07, 0.08)
  m <- 2

  i_eff <- (1 + j / m)^m - 1
  expected <- (1 + i_eff)^(-t)

  out <- yield_curve(
    t = t,
    i = j,
    i_type = "nominal_interest",
    m = m
  )

  expect_equal(out$i_effective[[1]], i_eff, tolerance = 1e-12)
  expect_equal(out$v[[1]], expected, tolerance = 1e-12)
})


test_that("yield_curve handles multiple curves in tidy list-columns", {
  curves <- tibble::tibble(
    curve_id = c("base", "stress"),
    t = list(c(1, 2, 3), c(1, 3, 5)),
    i = list(c(0.03, 0.04, 0.05), c(0.05, 0.055, 0.06))
  )

  out <- yield_curve(curves)

  expect_equal(nrow(out), 2L)
  expect_equal(out$curve_id, c("base", "stress"))
  expect_equal(
    out$v[[1]],
    (1 + c(0.03, 0.04, 0.05))^(-c(1, 2, 3)),
    tolerance = 1e-12
  )
  expect_equal(
    out$v[[2]],
    (1 + c(0.05, 0.055, 0.06))^(-c(1, 3, 5)),
    tolerance = 1e-12
  )
})


test_that("yield_curve supports custom tidy column names", {
  curves <- tibble::tibble(
    id = "A",
    maturity = list(c(1, 2, 3)),
    spot = list(c(0.04, 0.045, 0.05))
  )

  out <- yield_curve(
    curves,
    col_t = "maturity",
    col_i = "spot",
    .out = "discount",
    .keep = "used"
  )

  expect_named(out, c("maturity", "spot", "discount", "i_effective"))
  expect_equal(
    out$discount[[1]],
    (1 + curves$spot[[1]])^(-curves$maturity[[1]]),
    tolerance = 1e-12
  )
})


test_that("yield_curve .keep none returns only computed columns", {
  out <- yield_curve(
    t = c(1, 2),
    i = c(0.04, 0.05),
    .out = "df",
    .keep = "none"
  )

  expect_named(out, c("df", "i_effective"))
  expect_equal(out$df[[1]], c(1.04^-1, 1.05^-2), tolerance = 1e-12)
})


test_that("yield_curve plot output is a ggplot object", {
  out <- yield_curve(
    t = c(1, 2, 3),
    i = c(0.04, 0.045, 0.05),
    plot = TRUE
  )

  expect_true("yield_curve_plot" %in% names(out))
  expect_s3_class(out$yield_curve_plot[[1]], "ggplot")
})


test_that("yield_curve propagates missing values by row", {
  curves <- tibble::tibble(
    id = c("complete", "missing"),
    t = list(c(1, 2), c(1, 2)),
    i = list(c(0.04, 0.05), c(0.04, NA_real_))
  )

  out <- yield_curve(curves, .na = "propagate")

  expect_false(anyNA(out$v[[1]]))
  expect_true(all(is.na(out$v[[2]])))
  expect_true(all(is.na(out$i_effective[[2]])))
})


test_that("yield_curve can drop rows with missing required inputs", {
  curves <- tibble::tibble(
    id = c("complete", "missing"),
    t = list(c(1, 2), c(1, 2)),
    i = list(c(0.04, 0.05), c(0.04, NA_real_))
  )

  expect_warning(
    out <- yield_curve(curves, .na = "drop"),
    "Dropping row"
  )

  expect_equal(nrow(out), 1L)
  expect_equal(out$id, "complete")
})


test_that("yield_curve enforces a valid discrete maturity grid", {
  expect_error(
    yield_curve(t = c(1, 3, 2), i = c(0.04, 0.05, 0.055)),
    "strictly increasing"
  )

  expect_error(
    yield_curve(t = c(1, 2, 2), i = c(0.04, 0.05, 0.055)),
    "strictly increasing"
  )

  expect_error(
    yield_curve(t = c(1, 2), i = c(0.04, 0.05, 0.06)),
    "same length"
  )
})
