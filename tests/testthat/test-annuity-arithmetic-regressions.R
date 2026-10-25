test_that("arithmetic patterns retain their original financial basis", {
  for (force in c(-800, 800)) {
    for (pattern in c("increasing", "decreasing", "custom")) {
      payments <- switch(pattern, increasing = 1:4, decreasing = 4:1,
        custom = c(-3, -1, 1, 3))
      for (timing in c("immediate", "due")) {
        times <- (seq_len(4) - as.integer(timing == "due")) / 1600
        out <- annuity_arith(4 / 1600, k = 1600, i = force, i_type = "force",
          pattern = pattern, P1 = -3, g = 2, timing = timing, tidy = TRUE)
        expect_equal(out$present_value_factor, sum(payments * exp(-force * times)),
          tolerance = 1e-12)
        expect_equal(out$accumulated_value_factor,
          sum(payments * exp(force * (4 / 1600 - times))), tolerance = 1e-12)
      }
    }
  }
})

test_that("arithmetic monetary products recover before underflowing factors", {
  expect_equal(annuity_arith(1, i = 550 * log(2), i_type = "force", h = 1,
    pattern = "custom", P1 = 2^1000, g = 0) / 2^-100, 1, tolerance = 1e-12)
  out <- annuity_arith(2, i = -log(100), i_type = "force", h = 200,
    pattern = "custom", P1 = 0, g = 0, tidy = TRUE)
  expect_equal(out$present_value_factor, 0)
  expect_equal(out$accumulated_value_factor, 0)
})

test_that("arithmetic small rates and signed payments match explicit timelines", {
  for (force in c(-1e-12, 0, 1e-12, .04)) {
    times <- (1:12) / 4
    payments <- 10 - 3 * (0:11)
    out <- annuity_arith(3, k = 4, i = force, i_type = "force", h = 1,
      pattern = "custom", P1 = 10, g = -3, tidy = TRUE)
    expect_equal(out$present_value_factor, sum(payments * exp(-force * (1 + times))),
      tolerance = 1e-12)
    expect_equal(out$accumulated_value_factor, sum(payments * exp(force * (3 - times))),
      tolerance = 1e-12)
  }
  expect_true(is.na(annuity_arith(NA_real_, i = .04)))
  expect_true(is.na(annuity_arith(2, i = .04, i_type = NA_character_)))
})

test_that("arithmetic tabular scenarios preserve pattern-specific payment units", {
  tab <- tibble::tibble(term = c(2, 3), pattern = c("increasing", "custom"),
    first = c(NA_real_, 20), increment = c(NA_real_, -3))
  out <- tab |> dplyr::mutate(value = annuity_arith(term, i = 0,
    pattern = pattern, P1 = first, g = increment))
  expect_s3_class(out, "tbl_df")
  expect_equal(out$value, c(3, 51))
  expect_error(annuity_arith(2, i = .04, m = 0), "`m`")
  expect_error(annuity_arith(2, i = .04, pattern = "custom", g = NA_real_),
    "must not be NA")
})
