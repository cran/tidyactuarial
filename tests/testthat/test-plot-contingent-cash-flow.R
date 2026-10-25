contingent_values <- function(p, measure) {
  rows <- startsWith(as.character(p$data$measure), measure)
  p$data$value[rows]
}

test_that("survival plots retain independently valued payment rows", {
  lt <- tibble::tibble(x = 60:64, lx = c(1000, 800, 600, 300, 0))
  out <- apv_life_flow(lt, ages = 60, t = 0:3, cf = c(-100, 200, 300, 400), i = .04)
  p <- out |> plot_contingent_cash_flow(currency = "$")
  probabilities <- c(1, .8, .6, .3)
  expected <- c(-100, 200, 300, 400) * probabilities
  pv <- expected / 1.04^(0:3)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$data, "tbl_df")
  expect_equal(contingent_values(p, "Contractual"), out$cf)
  expect_equal(contingent_values(p, "Unconditional"), probabilities)
  expect_equal(contingent_values(p, "Expected cash"), expected)
  expect_equal(contingent_values(p, "APV"), pv)
  expect_identical(contingent_values(p, "APV"), out$pv)
  expect_equal(nrow(p$data), 16L)
  expect_identical((p + ggplot2::labs(title = "Custom"))$labels$title, "Custom")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("death intervals use unconditional probability and end of interval discounts", {
  # l60=1000,l61=800,l62=600,l63=300. Annual q61=.25, but death probability=.20.
  tab <- tibble::tibble(t = 1:3, cf = rep(10000, 3),
    probability = c(.2, .2, .3), discount = 1.04^-(1:3), event = "death_interval")
  p <- plot_contingent_cash_flow(tab)
  oracle <- 10000 * c(.2, .2, .3) / 1.04^(1:3)
  expect_equal(contingent_values(p, "Unconditional"), c(.2, .2, .3))
  expect_equal(contingent_values(p, "Expected cash"), c(2000, 2000, 3000))
  expect_equal(contingent_values(p, "APV"), oracle)
  lt <- data.frame(x = 60:64, lx = c(1000, 800, 600, 300, 0))
  expect_equal(sum(oracle), insurance_x(lt, x = 60, i = .04, n = 3,
    type = "term", benefit = 10000))
  expect_false(isTRUE(all.equal(sum(oracle), sum(10000 * c(.2, .25, .5) / 1.04^(1:3)))))
  expect_equal(nrow(plot_contingent_cash_flow(tab, view = "pv")$data), 3L)
})

test_that("annuity values and independent multi-life states agree", {
  lt <- data.frame(x = 60:64, lx = c(1000, 800, 600, 300, 0))
  out <- apv_life_flow(lt, 60, t = 1:3, cf = rep(1, 3), i = .04)
  expect_equal(sum(contingent_values(plot_contingent_cash_flow(out), "APV")),
    annuity_x(lt, x = 60, i = .04, n = 3, timing = "immediate"))
  for (status in c("first", "last", "reversionary")) {
    out <- apv_life_flow(lt, c(60, 61), t = 1:2, cf = c(100, 100),
      i = 0, status = status, alpha = .4)
    a <- c(.8, .6); b <- c(.75, .375)
    oracle <- switch(status, first = a*b, last = a+b-a*b,
      reversionary = a*b + .4*(a+b-2*a*b))
    p <- plot_contingent_cash_flow(out)
    expect_equal(contingent_values(p, "Expected cash"), 100 * oracle)
    if (status == "reversionary") {
      expect_true("Expected benefit fraction" %in% levels(p$data$measure))
      expect_false(any(grepl("probability", levels(p$data$measure))))
    } else expect_equal(contingent_values(p, "Unconditional"), oracle)
  }
  attr(out, "contingency") <- NULL
  expect_true("Payment weight" %in% levels(plot_contingent_cash_flow(out)$data$measure))
})

test_that("dates retain original actuarial times, simultaneous events remain separate", {
  lt <- data.frame(x = 60:64, lx = c(1000, 800, 600, 300, 0))
  dates <- as.Date(c("2026-07-01", "2027-01-01"))
  out <- apv_life_flow(lt, 60, date = dates, date0 = as.Date("2026-01-01"),
    cf = c(100, 200), i = .04)
  p <- plot_contingent_cash_flow(out, view = "pv")
  expect_identical(p$data$x, dates)
  expect_equal(p$data$t, as.numeric(dates - as.Date("2026-01-01")) / 365.25)
  expect_identical(p$data$value, out$pv)
  expect_no_warning(ggplot2::ggplot_build(p))
  tab <- tibble::tibble(t = c(1,1,1), cf = c(100,200,-50),
    probability = c(.8,.2,1), discount = rep(1/1.04,3),
    event = c("survival","death_interval","certain"))
  p <- plot_contingent_cash_flow(tab)
  expect_equal(contingent_values(p,"Unconditional"), c(.8,.2,1))
  expect_equal(contingent_values(p,"APV"), c(80,40,-50)/1.04)
  expect_equal(nrow(p$data), 12L)
})

test_that("stored rows are reused and malformed monetary identities are rejected", {
  tab <- tibble::tibble(t = 1, cf = 100, probability = .8,
    discount = .9, event = "survival", expected_cf = 80, pv = 72)
  attr(tab,"apv") <- -999
  p <- plot_contingent_cash_flow(tab)
  expect_match(p$labels$subtitle, "72.00", fixed = TRUE)
  tab$expected_cf <- 81
  expect_error(plot_contingent_cash_flow(tab), "expected_cf.*inconsistent")
  tab$expected_cf <- 80; tab$pv <- 73
  expect_error(plot_contingent_cash_flow(tab), "pv.*inconsistent")
  tab$cf <- 1e-250; tab$expected_cf <- 1e-250; tab$pv <- 1e-250
  expect_error(plot_contingent_cash_flow(tab), "inconsistent")
})

test_that("extreme products and zero events are finite without losing signs", {
  tab <- tibble::tibble(t = 1:3, cf = c(1e308, -1e308, 0),
    probability = c(1e-300, 0, 1), discount = c(1e-100, 1e300, 1e300),
    event = "survival")
  p <- plot_contingent_cash_flow(tab)
  expect_equal(contingent_values(p,"APV") / 1e-92, c(1,0,0), tolerance = 1e-10)
  expect_true(all(is.finite(p$data$plot_value)))
  expect_no_warning(ggplot2::ggplot_build(p))
  tab$cf <- c(1e-310,-1e-310,0); tab$probability <- c(1,1,0)
  tab$discount <- rep(1,3)
  p <- plot_contingent_cash_flow(tab)
  expect_true(all(is.finite(p$data$plot_value)))
  expect_true(all(contingent_values(p,"Contractual")[1:2] != 0))
})

test_that("invalid events, missing rows and summary objects fail explicitly", {
  tab <- tibble::tibble(t = 1, cf = 100, probability = .8, discount = .9, event = "survival")
  for (column in c("t","cf","probability","discount")) {
    bad <- tab; bad[[column]] <- NA_real_
    expect_error(plot_contingent_cash_flow(bad), "finite")
  }
  for (prob in c(-.1,1.1)) {
    bad <- tab; bad$probability <- prob
    expect_error(plot_contingent_cash_flow(bad), "probabilities")
  }
  bad <- tab; bad$event <- "certain"
  expect_error(plot_contingent_cash_flow(bad), "certain")
  bad$event <- "death_exact"
  expect_error(plot_contingent_cash_flow(bad), "event")
  expect_error(plot_contingent_cash_flow(tab[,1:4]), "event")
  expect_error(plot_contingent_cash_flow(tibble::tibble(apv = 100)), "payment-level")
  expect_error(plot_contingent_cash_flow(100), "payment-level")
  expect_error(plot_contingent_cash_flow(tab[FALSE,]), "nonempty")
  bad <- tab; bad$date <- as.Date(Inf, origin="1970-01-01")
  expect_error(plot_contingent_cash_flow(bad), "finite Date")
  expect_error(plot_contingent_cash_flow(tab, show_labels = NA), "show_labels")
  expect_error(plot_contingent_cash_flow(tab, label_digits = 51), "label_digits")
  expect_error(plot_contingent_cash_flow(tab, label_size = 0), "label_size")
})
