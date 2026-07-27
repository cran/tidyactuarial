test_that("plot_cash_flow returns a ggplot object", {
  p <- plot_cash_flow(
    C = c(-1000, 300, 400, 500),
    t = c(0, 1, 2, 3)
  )

  expect_s3_class(p, "ggplot")
})

test_that("positive-only cash flows do not reserve a negative region", {
  p <- plot_cash_flow(
    C = c(1000000, 1500000, 2000000),
    t = c(1, 2, 3)
  )

  built <- ggplot2::ggplot_build(p)
  y_range <- built$layout$panel_params[[1]]$y.range

  expect_gte(y_range[[1]], -1e-10)
  expect_gt(y_range[[2]], 0)
})

test_that("negative-only cash flows do not reserve a positive region", {
  p <- plot_cash_flow(
    C = c(-1000000, -1500000, -2000000),
    t = c(1, 2, 3)
  )

  built <- ggplot2::ggplot_build(p)
  y_range <- built$layout$panel_params[[1]]$y.range

  expect_lt(y_range[[1]], 0)
  expect_lte(y_range[[2]], 1e-10)
})

test_that("mixed cash-flow limits are proportional rather than forced symmetric", {
  p <- plot_cash_flow(
    C = c(-5000000, 1000000, 1500000, 2000000),
    t = c(0, 1, 2, 3)
  )

  built <- ggplot2::ggplot_build(p)
  y_range <- built$layout$panel_params[[1]]$y.range

  expect_gt(abs(y_range[[1]]), y_range[[2]])
})

test_that("zero cash flows produce a valid nondegenerate panel", {
  p <- plot_cash_flow(
    C = c(0, 0, 0),
    t = c(0, 1, 2)
  )

  built <- ggplot2::ggplot_build(p)
  y_range <- built$layout$panel_params[[1]]$y.range

  expect_lt(y_range[[1]], 0)
  expect_gt(y_range[[2]], 0)
})

test_that("large values use compact labels automatically", {
  p <- plot_cash_flow(
    C = c(-250000000, 50000000, 400000000),
    t = c(0, 1, 2),
    currency = "$"
  )

  expect_true(any(grepl("M|B|T", p$data$label)))
})

test_that("full labels remain available", {
  p <- plot_cash_flow(
    C = c(-250000000, 400000000),
    t = c(0, 1),
    label_format = "full"
  )

  expect_false(any(grepl("K|M|B|T", p$data$label)))
})

test_that("labels at horizontal extremes are aligned into the panel", {
  p <- plot_cash_flow(
    C = c(-1000000, 500000, 2500000),
    t = c(0, 1, 2)
  )

  expect_equal(p$data$label_hjust[[1]], 0)
  expect_equal(p$data$label_hjust[[3]], 1)
})

test_that("date inputs remain supported", {
  p <- plot_cash_flow(
    C = c(-1000, 450, 700),
    date = as.Date(c("2026-01-01", "2026-07-01", "2027-01-01")),
    i = 0.08
  )

  expect_s3_class(p, "ggplot")
})

test_that("deprecated argument names remain compatible", {
  p <- plot_cash_flow(
    payment = c(-1000, 300, 400, 500),
    time = c(0, 1, 2, 3),
    rate = 0.08
  )

  expect_s3_class(p, "ggplot")
})

convexity_cash_flow(cf = c(100,100,200,300),
                    t = 1:4,
                    rate = 0.06,
                    output = "value")
