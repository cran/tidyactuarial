test_that("path insurance uses original rates and monetary range", {
  d <- data.frame(Kx = c(0, 1, NA), Tx = c(.75, 1.5, NA))
  for (force in c(720, -1000)) for (timing in c("end_of_year", "moment_of_death")) {
    amount <- if (force > 0) 1e300 else 1e-300
    out <- mc_insurance(d, i = force, i_type = "force", benefit = amount,
                        type = "term", n = 1, timing = timing)
    date <- if (timing == "end_of_year") 1 else .75
    expect_equal(out$pv_benefit[1] / exp(log(amount) - force * date), 1,
                 tolerance = 1e-12)
    expect_identical(out$pv_benefit[2], 0)
    expect_true(is.na(out$pv_benefit[3]))
  }
})

test_that("path annuities discount monetary payments before losing range", {
  d <- data.frame(Kx = 1, Tx = 1.75)
  for (force in c(720, -1000)) for (timing in c("due", "immediate")) {
    amount <- if (force > 0) 1e300 else 1e-300
    out <- mc_annuity(d, i = force, i_type = "force", payment = amount,
                      type = "temporary", n = .5, k = 4, timing = timing)
    dates <- if (timing == "due") c(0, .25) else c(.25, .5)
    expected <- sum(exp(log(amount) - force * dates))
    expect_equal(out$pv_annuity / expected, 1, tolerance = 1e-12)
    expect_equal(out$n_payments, 2L)
  }
})

test_that("zero amounts and time zero retain their monetary meaning", {
  d <- data.frame(Kx = c(0, 1, NA), Tx = c(0, 1.5, NA))
  for (fun in c("mc_annuity", "mc_insurance")) {
    args <- list(.data = d, i = -1000, i_type = "force")
    args[[if (fun == "mc_annuity") "payment" else "benefit"]] <- 0
    out <- do.call(get(fun), args)
    value <- out[[if (fun == "mc_annuity") "pv_annuity" else "pv_benefit"]]
    expect_equal(value[1:2], c(0, 0))
    expect_true(is.na(value[3]))
  }
  # Log growth itself can overflow with an extreme valid nominal rate.
  args <- list(.data = d[1, ], i = -1e308 * (1 - 1e-15),
               i_type = "nominal_interest", m = 1e308)
  a <- do.call(mc_annuity, c(args, list(payment = 7, timing = "due")))
  z <- do.call(mc_insurance, c(args, list(benefit = 7, timing = "moment_of_death")))
  expect_equal(a$pv_annuity, 7)
  expect_equal(z$pv_benefit, 7)
})

test_that("all annuity types agree with independent dated cash flows", {
  d <- data.frame(Tx = c(0, .125, .5, .75, 1, 1.125, 1.5, 2, 2.75, 3.2, NA))
  d$Kx <- floor(d$Tx)
  for (type in c("whole", "temporary", "deferred", "deferred_temporary", "certain", "guaranteed"))
    for (k in c(1, 4)) for (timing in c("due", "immediate")) for (rate in c(-.02, 0, .05)) {
      h <- if (k == 1) 1 else .125
      args <- list(.data = d, i = rate, payment = 7, type = type, k = k, timing = timing)
      if (type %in% c("temporary", "deferred_temporary", "certain")) args$n <- 2
      if (type %in% c("deferred", "deferred_temporary")) args$h <- h
      if (type == "guaranteed") args$n_guar <- 2
      out <- do.call(mc_annuity, args)
      origin <- if (type %in% c("deferred", "deferred_temporary")) h else 0
      grid <- origin + (if (timing == "due") 0:20 else 1:20) / k
      expected <- lapply(seq_len(nrow(d)), function(j) {
        if (is.na(d$Tx[j]) && type != "certain") return(c(NA, NA, NA, NA))
        alive <- if (type == "certain") rep(TRUE, length(grid)) else if (k == 1)
          grid <= d$Kx[j] else if (timing == "due") grid < d$Tx[j] else grid <= d$Tx[j]
        if (type == "guaranteed") alive <- alive |
          if (timing == "due") grid < 2 else grid <= 2
        if (type %in% c("temporary", "deferred_temporary", "certain"))
          alive <- alive & if (timing == "due") grid < origin + 2 else grid <= origin + 2
        dates <- grid[alive]
        c(sum(7 * exp(-log1p(rate) * dates)), length(dates),
          if (length(dates)) min(dates) else NA,
          if (length(dates)) max(dates) else NA)
      })
      expected <- do.call(rbind, expected)
      expect_equal(out$pv_annuity, expected[, 1], tolerance = 1e-12)
      expect_equal(out$n_payments, expected[, 2])
      expect_equal(out$first_payment_time, expected[, 3])
      expect_equal(out$last_payment_time, expected[, 4])
    }
})

test_that("all insurance types agree with independent event ledgers", {
  d <- data.frame(Tx = c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5, NA))
  d$Kx <- floor(d$Tx)
  for (type in c("whole", "term", "deferred", "deferred_term", "pure_endowment", "endowment"))
    for (timing in c("end_of_year", "moment_of_death")) for (rate in c(-.02, 0, .05)) {
      args <- list(.data = d, i = rate, benefit = 11, type = type, timing = timing)
      if (type %in% c("term", "deferred_term", "pure_endowment", "endowment")) args$n <- 2
      if (type %in% c("deferred", "deferred_term")) args$h <- 1
      out <- do.call(mc_insurance, args)
      death <- if (timing == "end_of_year") d$Kx + 1 else d$Tx
      paid <- switch(type, whole = !is.na(death), term = death <= 2,
        deferred = death > 1, deferred_term = death > 1 & death <= 3,
        pure_endowment = if (timing == "end_of_year") d$Kx >= 2 else d$Tx > 2,
        endowment = !is.na(death))
      paid[is.na(death)] <- NA
      date <- if (type == "pure_endowment") rep(2, length(death)) else
        if (type == "endowment") pmin(death, 2) else death
      date <- ifelse(paid, date, NA_real_)
      expected <- ifelse(paid, 11 * exp(-log1p(rate) * date), 0)
      expect_equal(out$benefit_indicator, as.double(paid))
      expect_equal(out$benefit_time, date)
      expect_equal(out$pv_benefit, expected, tolerance = 1e-12)
    }
})

for (name in c("mc_annuity", "mc_insurance")) {
  test_that(paste(name, "retains large nominal conversion frequencies"), {
    out <- do.call(get(name), list(.data = data.frame(Kx = 1, Tx = 1.5),
      i = .05, i_type = "nominal_interest", m = 3e9))
    expect_equal(out$m, 3e9)
    expect_equal(out$i_effective, expm1(3e9 * log1p(.05 / 3e9)))
  })
  test_that(paste(name, "rejects conversion frequencies rounded to zero"), {
    expect_error(do.call(get(name), list(.data = data.frame(Kx = 1), i = .05, m = 1e-12)),
                 "m|positive|integer")
  })
}
test_that("annuity payment frequencies fail before integer overflow", {
  expect_error(mc_annuity(data.frame(Kx = 1, Tx = 1.5), i = .05, k = 3e9), "k|range")
  expect_error(mc_annuity(data.frame(Kx = 1, Tx = 1.5), i = .05, k = 1e-12), "k|positive|integer")
})

test_that("path valuation is deterministic and leaves the RNG state unchanged", {
  set.seed(81)
  before <- .Random.seed
  d <- data.frame(Kx = c(0, 1), Tx = c(.5, 1.5))
  a <- mc_annuity(d, i = .05, k = 4)
  z <- mc_insurance(d, i = .05)
  expect_identical(.Random.seed, before)
  expect_identical(mc_annuity(d, i = .05, k = 4), a)
  expect_identical(mc_insurance(d, i = .05), z)
})

test_that("C18 makes term and pure endowment a disjoint pathwise partition", {
  d <- data.frame(Tx = c(0, 1, 2-1e-12, 2, 2+1e-12, 3, NA))
  d$Kx <- floor(d$Tx)
  for (timing in c("end_of_year", "moment_of_death")) for (rate in c(-.02, 0, .05)) {
    values <- lapply(c("term", "pure_endowment", "endowment"), function(type)
      mc_insurance(d, i = rate, benefit = 10, type = type, n = 2, timing = timing))
    expect_equal(values[[1]]$pv_benefit + values[[2]]$pv_benefit,
                 values[[3]]$pv_benefit, tolerance = 1e-12)
    expect_equal(values[[1]]$benefit_indicator + values[[2]]$benefit_indicator,
                 values[[3]]$benefit_indicator)
  }
  pure <- mc_insurance(d, i = 0, benefit = 10, type = "pure_endowment",
                       n = 2, timing = "moment_of_death")
  expect_equal(pure$benefit_indicator, c(0, 0, 0, 0, 1, 1, NA))
  expect_equal(pure$pv_benefit, c(0, 0, 0, 0, 10, 10, NA))
  expect_equal(pure$benefit_time, c(NA, NA, NA, NA, 2, 2, NA))
})
