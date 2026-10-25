test_that("extreme diagrams retain finite geometry and readable nonzero labels", {
  p <- plot_cash_flow(C = c(-1.7e308, 1.7e308), t = 0:1)
  expect_true(all(is.finite(p$data$label_y)))
  expect_true(all(is.finite(p$coordinates$limits$y)))
  expect_true(all(nchar(p$data$label) < 25))
  tiny <- plot_cash_flow(C = 1e-20, t = 1)
  expect_match(tiny$data$label, "e-20", fixed = TRUE)
  expect_equal(p$data$C, c(-1.7e308, 1.7e308))
})

test_that("diagram inputs reject infinite dates and avoid integer coercion", {
  expect_error(plot_cash_flow(C = 1, date = as.Date(Inf, origin = "1970-01-01")),
               "finite.*date|date.*finite")
  expect_no_warning(plot_cash_flow(C = 100, t = 1, i = .05, m = 2^31))
})

test_that("dated PV follows the chosen day count and simultaneous signed flows", {
  dates <- as.Date(c("2028-01-01", "2029-01-01", "2029-01-01"))
  tab <- tibble::tibble(day = dates, cash = c(-100, 150, -40))
  for (denom in c(365, 360)) {
    p <- tab |> plot_cash_flow(C = cash, date = day, i = .04,
      day_count = paste0("act/", denom), label_digits = 8)
    pv <- -100 + 110 * 1.04^(-366 / denom)
    expect_match(p$labels$subtitle, formatC(pv, format = "f", digits = 8), fixed = TRUE)
    expect_equal(p$data$C, c(-100, 110))
    expect_s3_class(p$data$x, "Date")
  }
  separate <- tab |> plot_cash_flow(C = cash, date = day, aggregate = FALSE)
  expect_equal(nrow(separate$data), 3L)
  expect_equal(separate$data$x, dates)
})

test_that("dense plots preserve all flows and allow labels to be disabled", {
  tab <- tibble::tibble(C = rep(c(-10, 20), 60), t = seq_len(120))
  p <- tab |> plot_cash_flow(C = C, t = t)
  expect_equal(nrow(p$data), 120L)
  expect_true(p$layers[[3]]$geom_params$check_overlap)
  expect_length(plot_cash_flow(tab, C = C, t = t, show_labels = FALSE)$layers, 2L)
  expect_s3_class(p + ggplot2::theme_bw(), "ggplot")
})
