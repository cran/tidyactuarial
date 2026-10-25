test_that("C06 propagates unknown portfolio identifiers without silently dropping rows", {
  for (kind in c("duration", "convexity")) {
    fn <- get(paste0("portfolio_", kind), asNamespace("tidyactuarial"))
    metric <- if (kind == "duration") "D" else "C"
    output <- paste0(metric, "_P")
    for (ids in list(c("b", NA_character_, "a", NA_character_),
                     c(2, NA_real_, 1, NA_real_),
                     factor(c("b", NA_character_, "a", NA_character_),
                            levels = c("a", "b", "unused")))) {
      args <- list(P = c(100, 200, 300, 400), portfolio_id = ids)
      args[[metric]] <- c(2, 4, 6, 8)
      out <- do.call(fn, args)
      expect_equal(nrow(out), 3L)
      expect_identical(class(out$portfolio_id), class(ids))
      expect_identical(is.na(out$portfolio_id), c(FALSE, FALSE, TRUE))
      expect_equal(out[[output]][1:2], c(6, 2))
      expect_equal(out$P_total[1:2], c(300, 100))
      expect_true(is.na(out[[output]][3]))
      expect_true(is.na(out$P_total[3]))
      expect_true(is.na(out$n_positions[3]))
      expect_error(do.call(fn, c(args, list(.na = "error"))), "Missing values")
      expect_warning(dropped <- do.call(fn, c(args, list(.na = "drop"))), "Dropping")
      expect_equal(nrow(dropped), 2L)
      expect_false(anyNA(dropped$portfolio_id))
    }
  }
})

test_that("C06 handles all-missing identifiers and a literal NA label in tidy input", {
  for (kind in c("duration", "convexity")) {
    fn <- get(paste0("portfolio_", kind), asNamespace("tidyactuarial"))
    metric <- if (kind == "duration") "D" else "C"
    data <- tibble::tibble(book = c("NA", NA_character_),
                          value = c(100, 200), sensitivity = c(3, 7))
    args <- list(.data = data, col_portfolio = "book", col_P = "value",
                 .out = "risk", .out_value = "value_total", .out_n = "count")
    args[[paste0("col_", metric)]] <- "sensitivity"
    out <- do.call(fn, args)
    expect_named(out, c("book", "risk", "value_total", "count"))
    expect_identical(out$book, c("NA", NA_character_))
    expect_equal(out$risk[1], 3)
    expect_true(is.na(out$risk[2]))
    data$book[] <- NA_character_
    args$.data <- data
    out <- do.call(fn, args)
    expect_equal(nrow(out), 1L)
    expect_true(is.na(out$book))
    expect_true(all(is.na(out$risk), is.na(out$value_total), is.na(out$count)))
    expect_warning(empty <- do.call(fn, c(args, list(.na = "drop"))), "Dropping")
    expect_equal(nrow(empty), 0L)
  }
})
