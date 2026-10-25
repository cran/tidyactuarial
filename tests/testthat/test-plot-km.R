make_km_plot_object <- function() {
  km_lifetable(
    time = c(1, 2, 3, 4, 5, 6),
    status = c(1, 0, 1, 1, 0, 1),
    breaks = 0:6
  )
}

test_that("plot_km accepts the full km_lifetable output", {
  out <- make_km_plot_object()

  p <- plot_km(out)

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 2L)
})

test_that("plot_km can omit the confidence ribbon", {
  out <- make_km_plot_object()

  p <- plot_km(
    out$km,
    conf_int = FALSE,
    title = "Empirical survival"
  )

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 1L)
  expect_equal(p$labels$title, "Empirical survival")
})

test_that("plot_km supports custom time, survival and confidence columns", {
  custom <- tibble::tibble(
    tt = c(1, 2, 3),
    ss = c(0.9, 0.8, 0.7),
    lo = c(0.8, 0.7, 0.6),
    hi = c(0.98, 0.9, 0.82)
  )

  p <- plot_km(
    custom,
    time_col = "tt",
    surv_col = "ss",
    lower_col = "lo",
    upper_col = "hi",
    conf_int = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 2L)
})

test_that("plot_km gracefully plots without a ribbon when CI columns are absent", {
  km <- tibble::tibble(
    time = c(1, 2, 3),
    S = c(0.9, 0.8, 0.7)
  )

  p <- plot_km(km, conf_int = TRUE)

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 1L)
})

test_that("plot_km validates the required plotting columns", {
  expect_error(
    plot_km(
      tibble::tibble(time = 1:3, S = c(0.9, 0.8, 0.7)),
      time_col = "missing_time"
    ),
    "not found"
  )

  expect_error(
    plot_km(
      tibble::tibble(time = 1:3, S = c(0.9, 0.8, 0.7)),
      surv_col = "missing_survival"
    ),
    "not found"
  )
})
