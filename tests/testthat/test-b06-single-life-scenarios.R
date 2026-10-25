test_that("annuity scenarios recycle scalars and retain independent row meanings", {
  lt <- data.frame(x = 60:65, lx = c(1000, 800, 500, 200, 50, 0))
  ages <- c(60, 61, 60)
  rate <- c(.05, .03, .06)
  basis <- c("effective", "force", "nominal_interest")
  m <- c(1, 1, 12)
  n <- c(2, 1.5, 1)
  h <- c(0, 0, 1)
  k <- c(1, 4, 12)
  timing <- c("due", "immediate", "due")
  delta <- c(log(1.05), .03, 12 * log(1 + .06 / 12))
  expected <- vapply(1:3, function(z) {
    j <- if (timing[z] == "due") 0:(n[z] * k[z] - 1) else 1:(n[z] * k[z])
    t <- h[z] + j / k[z]
    a <- ages[z] + floor(t)
    s <- t - floor(t)
    l <- (1 - s) * lt$lx[match(a, lt$x)] + s * lt$lx[match(a + 1, lt$x)]
    sum(l / lt$lx[match(ages[z], lt$x)] * exp(-delta[z] * t)) / k[z]
  }, numeric(1))
  out <- annuity_x(lt, ages, rate, i_type = basis, m = m, n = n, h = h,
                   k = k, timing = timing, tidy = TRUE)
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 3)
  expect_equal(out$x, ages)
  expect_equal(out$apv, expected)
  expect_equal(annuity_x(lt, ages, rate, i_type = basis, m = m, n = n,
                         h = h, k = k, timing = timing), expected)
  for (j in 1:3) {
    scalar <- annuity_x(lt, ages[j], rate[j], i_type = basis[j], m = m[j],
                        n = n[j], h = h[j], k = k[j], timing = timing[j], tidy = TRUE)
    expect_equal(out[j, ], scalar)
  }
})

test_that("insurance scenarios support contract types, benefits and nominal bases", {
  lt <- data.frame(x = 60:65, lx = c(1000, 800, 500, 200, 50, 0))
  x <- c(60, 61, 60)
  n <- c(Inf, 2, 2)
  h <- c(0, 0, 1)
  type <- c("whole", "term", "endowment")
  benefit <- c(100, -200, 300)
  rate <- c(.05, .03, .06)
  basis <- c("effective", "force", "nominal_discount")
  m <- c(1, 1, 12)
  delta <- c(log(1.05), .03, -12 * log(1 - .06 / 12))
  expected <- vapply(1:3, function(j) {
    years <- if (type[j] == "whole") 66 - x[j] - h[j] else n[j]
    t <- h[j] + 1:years
    counts <- function(age) ifelse(age == 66, 0, lt$lx[match(age, lt$x)])
    deaths <- counts(x[j] + t - 1) - counts(x[j] + t)
    value <- sum(deaths * exp(-delta[j] * t))
    if (type[j] == "endowment") {
      value <- value + counts(x[j] + h[j] + years) * exp(-delta[j] * (h[j] + years))
    }
    benefit[j] * value / counts(x[j])
  }, numeric(1))
  out <- insurance_x(lt, x, rate, i_type = basis, m = m, n = n, h = h,
                     type = type, benefit = benefit, tidy = TRUE)
  expect_s3_class(out, "tbl_df")
  expect_equal(out$apv, expected)
  expect_equal(out$type, type)
  expect_equal(insurance_x(lt, x, rate, i_type = basis, m = m, n = n, h = h,
                           type = type, benefit = benefit), expected)
  for (j in 1:3) {
    scalar <- insurance_x(lt, x[j], rate[j], i_type = basis[j], m = m[j],
                          n = n[j], h = h[j], type = type[j], benefit = benefit[j], tidy = TRUE)
    expect_equal(out[j, ], scalar)
  }
})

test_that("single-life scenarios work in mutate and with inherited contracts", {
  lt <- data.frame(x = 60:65, lx = c(1000, 800, 500, 200, 50, 0))
  policies <- tibble::tibble(age = c(60, 61, 62), rate = c(0, .03, .05))
  out <- dplyr::mutate(policies,
    a = annuity_x(lt, x = age, i = rate, n = 2, timing = "due"),
    b = insurance_x(lt, x = age, i = rate, n = 2, type = "term"))
  expect_equal(out$a, c(1.8, 1 + .625 / 1.03, 1 + .4 / 1.05))
  expect_equal(out$b, c(.5, (.375 / 1.03 + .375 / 1.03^2),
                        .6 / 1.05 + .3 / 1.05^2))
  contract <- life_contract(lt, x = 60, i = .04, i_type = "force")
  expect_equal(contract |> annuity_x(n = c(1, 2)),
               annuity_x(lt, 60, .04, i_type = "force", n = c(1, 2)))
  expect_equal(contract |> insurance_x(n = c(1, 2), type = "term"),
               insurance_x(lt, 60, .04, i_type = "force", n = c(1, 2), type = "term"))
  # Defaults, including NULL whole-life n and frac, are shared scalar controls.
  expect_length(annuity_x(lt, x = c(60, 61), i = .05), 2)
  expect_length(insurance_x(lt, x = c(60, 61), i = .05), 2)
  expect_length(insurance_x(mortality_table = lt, age = c(60, 61), rate = .05,
                            term_years = 2, insurance_type = "term"), 2)
})

test_that("scenario sizes and invalid rows fail without silent recycling", {
  lt <- data.frame(x = 60:65, lx = c(1000, 800, 500, 200, 50, 0))
  for (fun in list(annuity_x, insurance_x)) {
    expect_error(fun(lt, x = c(60, 61), i = c(.01, .02, .03)), "common positive length")
    expect_error(fun(lt, x = numeric(), i = .05), "positive length")
    expect_error(fun(lt, x = c(60, NA), i = .05), "Scenario 2:.*x")
    expect_error(fun(lt, x = c(60, 61), i = c(.05, NA)), "Scenario 2:.*i")
    expect_error(fun(lt, x = c(60, 61), i = .05, tidy = c(TRUE, FALSE)), "logical scalar")
  }
})

test_that("C10 pays zero-term endowments at the deferred maturity", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 50, 0))
  expect_equal(insurance_x(lt, 60, .05, n = 0, type = "endowment", benefit = 100), 100)
  expect_equal(insurance_x(lt, 60, .05, n = 0, h = 1,
                           type = "endowment", benefit = 100), 80 / 1.05)
  expect_equal(insurance_x(lt, 60, .05, n = 0, type = "term", benefit = 100), 0)
  out <- insurance_x(lt, 60, .05, n = c(0, 1), type = "endowment", tidy = TRUE)
  expect_equal(out$apv, c(1, 1 / 1.05))
  expect_equal(out$pure_endowment_factor, c(1, .8 / 1.05))
})

test_that("C11 gives empty annuities the ordinary diagnostic schema", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 50, 0))
  empty <- annuity_x(lt, 60, .05, n = 0, h = 1, tidy = TRUE)
  ordinary <- annuity_x(lt, 60, .05, n = 1, h = 1, tidy = TRUE)
  expect_identical(names(empty), names(ordinary))
  expect_equal(empty$apv, 0)
  expect_equal(empty$annuity_value_at_start, 0)
  expect_equal(empty$pure_endowment_factor, 1)
  expect_equal(empty$deferment_factor, .8 / 1.05)
  out <- annuity_x(lt, 60, .05, n = c(0, 1), tidy = TRUE)
  expect_equal(out$annuity_value_at_start, c(0, .8 / 1.05))
  expect_equal(out$pure_endowment_factor, c(1, .8 / 1.05))
})
