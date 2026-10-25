test_that("cash-flow plots preserve tabular pipe and modifiable ggplot contracts", {
  tab <- tibble::tibble(time = c(0, 1, 1), amount = c(-100, 50, 60))
  bare <- tab |> plot_cash_flow(C = amount, t = time, i = .1)
  quoted <- tab |> plot_cash_flow(C = "amount", t = "time", i = .1)
  expect_s3_class(bare, "ggplot")
  expect_equal(bare$data$C, c(-100, 110))
  expect_equal(bare$data, quoted$data)
  expect_identical((bare + ggplot2::labs(title = "Custom"))$labels$title, "Custom")
})

test_that("plot present values retain the original force and representable products", {
  tab <- tibble::tibble(time = 1, amount = 2^1000)
  out <- tab |> plot_cash_flow(C = amount, t = time, i = 1100 * log(2),
    i_type = "force", normalize = TRUE, label_format = "full", label_digits = 40)
  expected <- 2^-100
  expected_text <- formatC(expected, format = "f", digits = 40)
  expect_true(grepl(paste0("PV = ", expected_text), out$labels$subtitle, fixed = TRUE))
  expect_equal(out$data$height, 1)
})

test_that("plotting preserves negative timeline valuation and dated aggregation", {
  out <- plot_cash_flow(C = 100, t = -1, i = .1)
  expect_true(grepl("PV = 110.00", out$labels$subtitle, fixed = TRUE))
  tab <- tibble::tibble(day = as.Date(c("2026-01-01", "2027-01-01", "2027-01-01")),
    cash = c(-100, 50, 60))
  dated <- tab |> plot_cash_flow(C = cash, date = day, i = .1)
  expect_s3_class(dated$data$x, "Date")
  expect_equal(dated$data$C, c(-100, 110))
})

test_that("unrepresentable graph summaries are rejected explicitly", {
  expect_error(plot_cash_flow(C = c(1e308, 1e308), t = c(0, 0)),
    "Aggregated cash flows")
  expect_error(plot_cash_flow(C = 1e308, t = 1, i = -1, i_type = "force"),
    "Computed present value")
})
