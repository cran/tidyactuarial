test_that("contract storage does not turn finite integer-valued inputs into NA", {
  lt <- tibble::tibble(x = 40:46, lx = 1000 * 0.8^(0:6))
  big <- 2^31
  for (field in c("x", "m")) {
    args <- list(lt = lt, x = 40, i = 0.03)
    args[[field]] <- big
    contract <- do.call(life_contract, args)
    expect_identical(contract[[field]], big)
    expect_identical(contract$age, contract$x)
  }
  contract <- life_contract(list(lt, lt), "joint", x = -big, y = big,
                            i = 0.03, m = big)
  expect_identical(contract$x, -big)
  expect_identical(contract$y, big)
  expect_identical(contract$age_x, contract$x)
  expect_identical(contract$age_y, contract$y)
  base <- life_contract(lt, x = 40, i = 0.03)
  for (field in c("n", "h", "k")) {
    args <- list(contract = base, type = "term", n = 3)
    args[[field]] <- big
    result <- do.call(add_insurance, args)
    expect_identical(result$insurance[[field]], big)
    expect_identical(result$lt, lt)
  }
  schedule <- add_premium_schedule(base, k = big, n_prem = 1 / big)
  expect_identical(schedule$premium_schedule$k, big)
  expect_identical(schedule$premium_schedule$n_prem, 1 / big)
  expect_null(base$insurance)
  expect_null(base$premium_schedule)
})

test_that("contract rates respect their financial domain without annualizing storage", {
  lt <- data.frame(x = 40:46, lx = 1000 * 0.8^(0:6))
  expect_error(life_contract(lt, x = 40, i = -1), "greater than -1")
  expect_error(life_contract(lt, x = 40, i = -12,
                             i_type = "nominal_interest", m = 12), "positive")
  expect_error(life_contract(lt, x = 40, i = 12,
                             i_type = "nominal_discount", m = 12), "positive")
  for (rate in c(-1000, 1000)) {
    contract <- life_contract(lt, x = 40, i = rate, i_type = "force")
    expect_identical(contract$i, rate)
    expect_identical(contract$rate, rate)
    expect_identical(contract$i_type, "force")
  }
  contract <- life_contract(lt, x = 40, i = 0.12,
                            i_type = "nominal_interest", m = 12)
  expect_identical(contract$i, 0.12)
  expect_identical(contract$rate, 0.12)
  expect_identical(contract$m, 12L)
  big <- 2^31
  for (type in c("nominal_interest", "nominal_discount")) {
    contract <- life_contract(lt, x = 40, i = 0.06, i_type = type, m = big)
    growth <- if (type == "nominal_interest")
      big * log1p(0.06 / big) else -big * log1p(-0.06 / big)
    expect_identical(contract$m, big)
    expect_identical(contract$rate_type, type)
    expect_equal(annuity_x(contract, n = 2, timing = "due"),
                 1 + 0.8 * exp(-growth))
  }
})

test_that("premium grids reject overflow and empty rounded schedules explicitly", {
  lt <- data.frame(x = 40:46, lx = 1000 * 0.8^(0:6))
  base <- life_contract(lt, x = 40, i = 0.03)
  expect_error(add_premium_schedule(base, k = 2, n_prem = 1e308),
               "finite.*positive.*integer")
  expect_error(add_premium_schedule(base, n_prem = 1e-12),
               "finite.*positive.*integer")
  expect_error(add_premium_schedule(base, k = 2, n_prem = 0.3),
               "integer")
  expect_error(add_insurance(base, type = "term", n = 1e-12), "integer")
  expect_error(add_insurance(base, type = "endowment", n = 0), "positive")
  expect_null(add_premium_schedule(base)$premium_schedule$n_prem)
  expect_identical(add_premium_schedule(base, n_prem = Inf)$premium_schedule$n_prem,
                   Inf)
  expect_identical(add_premium_schedule(base, k = 12, n_prem = 2.5)$premium_schedule$n_prem,
                   2.5)
})

test_that("stored components produce the independent premium cash ledger", {
  lt <- tibble::tibble(x = 40:46, lx = 1000 * 0.8^(0:6))
  base <- life_contract(lt, x = 40, i = 0.06,
                        i_type = "nominal_interest", m = 2)
  contract <- base |>
    add_insurance(type = "term", benefit = 500, n = 3, h = 1, frac = "UDD") |>
    add_premium_schedule(k = 2, n_prem = 2)
  growth <- 2 * log1p(0.06 / 2)
  benefit_apv <- sum(500 * 0.8^(1:3) * 0.2 * exp(-growth * (2:4)))
  times <- (0:3) / 2
  survival <- 0.8^floor(times) * (1 - (times - floor(times)) * 0.2)
  premium_factor <- sum(survival * exp(-growth * times)) / 2
  expect_equal(premium_x(contract, output = "value"), benefit_apv / premium_factor)
  expect_identical(contract$lt, lt)
  expect_identical(contract$mortality_table, lt)
  expect_identical(contract$insurance$n, 3L)
  expect_identical(contract$insurance$h, 1L)
  expect_identical(contract$insurance$k, 1L)
  expect_identical(contract$premium_schedule$k, 2L)
  expect_identical(base$insurance, NULL)
  # Explicit consumer inputs override the stored specification.
  expect_equal(premium_x(contract, benefit = 1000, output = "value"),
               2 * benefit_apv / premium_factor)
  expect_equal(premium_x(contract, i = 0, i_type = "effective", output = "value"),
               sum(500 * 0.8^(1:3) * 0.2) / (sum(survival) / 2))
})

test_that("component assembly preserves snapshot assumptions and other fields", {
  lt <- tibble::tibble(x = 40:46, lx = 1000 * 0.8^(0:6))
  base <- life_contract(lt, x = 40, i = 0.03)
  base$note <- "retain this user field"
  first <- base |> add_insurance(frac = "CML") |> add_premium_schedule()
  second <- base |> add_premium_schedule() |> add_insurance(frac = "CML")
  expect_identical(first$insurance$frac, "CF")
  expect_identical(first$premium_schedule$frac, "CF")
  expect_identical(second$premium_schedule$frac, "UDD")
  expect_identical(second$insurance$frac, "CF")
  expect_identical(first$note, base$note)
  expect_identical(second$note, base$note)
  expect_identical(first$insurance$n, Inf)
  expect_identical(names(first$insurance),
                   c("type", "benefit", "n", "h", "k", "frac", "status"))
  expect_identical(names(first$premium_schedule),
                   c("k", "timing", "premium_start", "n_prem", "woolhouse", "frac"))
  calls <- 0L
  benefit <- function(t) { calls <<- calls + 1L; rep(100, length(t)) }
  variable <- add_insurance(base, type = "variable_k", benefit = benefit,
                             n = 1.5, k = 2)
  expect_identical(calls, 0L)
  expect_identical(variable$insurance$benefit, benefit)
  expect_identical(variable$insurance$n, 1.5)
  expect_identical(variable$insurance$k, 2L)
  inferred <- add_insurance(base, type = "variable_k", benefit = c(100, 200))
  expect_identical(inferred$insurance$n, Inf)
})
