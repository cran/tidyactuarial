test_that("factors retain the original force when annual equivalents saturate", {
  cases <- list(
    list(type = "force", i = -40, m = 1, t = 1, a = exp(-40)),
    list(type = "force", i = 800, m = 1, t = .001, a = exp(.8)),
    list(type = "force", i = -800, m = 1, t = .001, a = exp(-.8)),
    # One conversion period: nominal rates imply these exact factors.
    list(type = "nominal_interest", i = -500, m = 1000, t = .001, a = .5),
    list(type = "nominal_interest", i = 2000, m = 1000, t = .001, a = 3),
    list(type = "nominal_discount", i = 500, m = 1000, t = .001, a = 2)
  )
  for (z in cases) {
    args <- list(i = z$i, i_type = z$type, m = z$m, t = z$t)
    expect_equal(do.call(accumulation_factor, args) / z$a, 1, tolerance = 1e-12)
    expect_equal(do.call(discount_factor, args) * z$a, 1, tolerance = 1e-12)
    expect_equal(do.call(discount_factor_spot, args) * z$a, 1, tolerance = 1e-12)
    expect_equal(do.call(future_value, c(list(C = 2), args)) / z$a, 2,
                 tolerance = 1e-12)
    expect_equal(do.call(present_value, c(list(C = 2), args)) * z$a, 2,
                 tolerance = 1e-12)
    expect_equal(do.call(pv_flow, c(list(cf = c(1, 2)),
      modifyList(args, list(t = c(0, z$t))))) / (1 + 2 / z$a), 1,
      tolerance = 1e-12)
    expect_equal(do.call(fv_flow, c(list(cf = c(1, 2)),
      modifyList(args, list(t = c(0, z$t))))) / (z$a + 2), 1,
      tolerance = 1e-12)
  }
})

test_that("rate-equivalence tables retain representable values from an extreme force", {
  out <- interest_equivalents(i_type = "force", i = -40, m = 1)
  get <- function(f) out$value[out$family == f]
  expect_equal(get("force"), -40)
  expect_equal(get("discount_factor") / exp(40), 1, tolerance = 1e-12)
  expect_equal(get("discount") / (1 - exp(40)), 1, tolerance = 1e-12)
  expect_equal(get("nominal_discount") / (1 - exp(40)), 1, tolerance = 1e-12)
})

test_that("cash-flow ratios avoid overflow of common accumulation levels", {
  # Discount and accumulation levels at times 1000 and 1001 overflow in
  # opposite directions, but one year of interest remains exp(1).
  expect_equal(fv_flow(cf = c(100, 50), t = c(1000, 1001),
    i = 1, i_type = "force"), 100 * exp(1) + 50)
  expect_equal(fv_flow(cf = c(100, 50), t = c(1000, 1001),
    i = -1, i_type = "force"), 100 / exp(1) + 50)
  expect_equal(fv_flow(cf = c(100, 50), t = c(1000, 1000),
    i = 1e308, i_type = "force"), 150)
  expect_error(fv_flow(cf = c(1, 1), t = c(.001, .001),
    i = c(800, 801), i_type = "force"), "spot rates must agree")
})

test_that("interval discounting does not require a representable accumulation", {
  # exp(740) overflows, while exp(-740) is a representable subnormal.
  expected <- exp(-740)
  expect_equal(discount_factor(t = 1, i = 740, i_type = "force") / expected, 1)
  expect_equal(discount_factor(t = 1, delta = function(t) 740) / expected, 1)
  a <- function(t) ifelse(t == 0, 1e-200, 1e110)
  expect_equal(discount_factor(t = 1, a = a) / 1e-310, 1, tolerance = 1e-10)
  expect_error(accumulation_factor(t = 1, i = 740, i_type = "force"),
    "positive and finite")
})

test_that("single-payment tiny-rate results agree with the nontrivial limit", {
  expect_equal(present_value(C = 2, i = 1e-20, t = 1e20), 2 / exp(1))
  expect_equal(future_value(C = 2, i = 1e-20, t = 1e20), 2 * exp(1))
  expect_equal(pv_flow(cf = 2, i = 1e-20, t = 1e20), 2 / exp(1))
  expect_equal(fv_flow(cf = c(2, 0), i = 1e-20, t = c(0, 1e20)), 2 * exp(1))
})
