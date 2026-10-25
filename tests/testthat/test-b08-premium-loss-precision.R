test_that("premium amounts survive loss of range in the annuity coefficient", {
  for (amount in c(1e-300, 1e300)) {
    basis <- 1 / amount
    d <- data.frame(pv_benefit = amount, pv_annuity = basis)
    p <- mc_premium(d, k = 1, annuity_payment = basis)
    expect_equal(p$premium_per_payment / amount, 1, tolerance = 1e-12)
    expect_equal(p$premium_annualized / amount, 1, tolerance = 1e-12)
    expect_lt(abs(p$mc_equivalence_residual) / amount, 1e-12)
    loss <- mc_loss(p, annuity_payment = basis)
    expect_equal(loss$pv_premiums / amount, 1, tolerance = 1e-12)
    expect_lt(abs(loss$L) / amount, 1e-12)
  }
})

test_that("all premium units retain representable pathwise income", {
  for (amount in c(1e-300, 1e300)) for (k in c(1, 4, 12)) {
    basis <- 1 / amount
    d <- data.frame(pv_benefit = amount * c(1, 2, NA),
                    pv_annuity = basis * c(1, 2, NA))
    for (unit in c("annualized", "per_payment")) {
      out <- mc_loss(d, P = if (unit == "annualized") amount * k else amount,
                     premium_unit = unit, k = k, annuity_payment = basis)
      expect_equal(out$pv_premiums / amount, c(1, 2, NA), tolerance = 1e-12)
      expect_equal(out$L[1:2] / amount, c(0, 0), tolerance = 1e-12)
      expect_true(is.na(out$L[3]))
      expect_equal(out$premium_per_payment / amount, rep(1, 3), tolerance = 1e-12)
    }
  }
})

test_that("annualization recovers a quota below double range", {
  d <- data.frame(pv_benefit = 1e-200, pv_annuity = 1e200)
  out <- mc_premium(d, k = 1e9, annuity_payment = 1e70)
  # Exact q=1e-400 and installment=1e-330; annual amount=1e-321.
  expect_equal(out$premium_per_payment, 0)
  expect_equal(out$premium_annualized / 1e-321, 1, tolerance = .01)
})

test_that("subnormal means do not erase positive paired exposure", {
  tiny <- .Machine$double.xmin * .Machine$double.eps
  d <- data.frame(pv_benefit = c(tiny, 0, 0), pv_annuity = c(tiny, 0, 0))
  out <- mc_premium(d)
  expect_equal(out$premium_annuity_scale, rep(1, 3))
  expect_equal(out$premium_per_payment, rep(1, 3))
  d$pv_benefit <- c(2 * tiny, 0, 0)
  d$pv_annuity <- c(4 * tiny, 0, 0)
  out <- mc_premium(d)
  expect_equal(out$premium_annuity_scale, rep(.5, 3), tolerance = 1e-12)
  expect_equal(out$premium_per_payment, rep(.5, 3), tolerance = 1e-12)
})

test_that("premium equivalence uses paired means and contract-specific groups", {
  d <- data.frame(g = rep(c("a", "b", NA), each = 4),
    pv_benefit = rep(c(2, 10, NA, 1000), 3),
    pv_annuity = rep(c(1, 3, 1000, NA), 3),
    payment_per_payment = rep(c(.25, 2, .5), each = 4),
    payments_per_year = rep(c(4, 2, 12), each = 4))
  out <- mc_premium(d, by = "g")
  expect_equal(out$premium_annuity_scale, rep(3, 12))
  expect_equal(out$premium_per_payment, 3 * d$payment_per_payment)
  expect_equal(out$premium_annualized, 3 * d$payment_per_payment * d$payments_per_year)
  expect_equal(out$mc_equivalence_residual, rep(0, 12))
  grouped <- mc_premium(dplyr::group_by(d, g))
  expect_equal(grouped, tibble::as_tibble(out))
  expect_false(dplyr::is_grouped_df(grouped))
  expect_error(mc_premium(d, by = "g", na_rm = FALSE), "Missing")
  expect_error(mc_premium(data.frame(pv_benefit = 0, pv_annuity = 0)), "positive")
})

test_that("explicit path books give equivalent premiums and losses across units", {
  rate <- .04
  benefit <- c(100, 0, 100, 0)
  death_date <- c(1, 2, 3, 4)
  Z <- benefit * (1 + rate)^(-death_date)
  for (k in c(1, 2, 4, 12)) {
    dates <- lapply(death_date, function(t) (0:(t * k - 1)) / k)
    factors <- vapply(dates, function(t) sum((1 + rate)^(-t)), 0)
    quota <- sum(Z) / sum(factors)
    for (c in c(1 / k, 1, 7)) {
      d <- data.frame(pv_benefit = Z, pv_annuity = c * factors,
                      payment_per_payment = c, payments_per_year = k)
      premium <- mc_premium(d)
      expect_equal(premium$premium_per_payment, rep(quota, 4), tolerance = 1e-12)
      expect_equal(premium$premium_annualized, rep(k * quota, 4), tolerance = 1e-12)
      for (unit in c("annualized", "per_payment", "annuity_scale")) {
        P <- switch(unit, annualized = quota * k, per_payment = quota, annuity_scale = quota / c)
        loss <- mc_loss(d, P = P, premium_unit = unit)
        expect_equal(loss$pv_premiums, quota * factors, tolerance = 1e-12)
        expect_equal(loss$L, Z - quota * factors, tolerance = 1e-12)
        expect_lt(abs(mean(loss$L)), 1e-12)
      }
    }
  }
})

for (name in c("mc_premium", "mc_loss")) {
  test_that(paste(name, "retains a large valid payment frequency"), {
    args <- list(.data = data.frame(pv_benefit = 2, pv_annuity = 1),
                 k = 3e9, annuity_payment = 1)
    if (name == "mc_loss") args$P <- 2
    out <- do.call(get(name), args)
    expect_equal(out$payments_per_year, 3e9)
    expect_equal(out$premium_annualized, 6e9)
  })
}

test_that("zero premiums and missing paths retain their meanings", {
  d <- data.frame(pv_benefit = c(0, 0, NA), pv_annuity = c(1e300, 0, NA))
  p <- mc_premium(d, annuity_payment = 1e300)
  expect_equal(p$premium_per_payment, c(0, 0, 0))
  for (unit in c("annualized", "per_payment", "annuity_scale")) {
    out <- mc_loss(d, P = 0, premium_unit = unit, annuity_payment = 1e-300)
    expect_equal(out$pv_premiums, c(0, 0, NA))
    expect_equal(out$L, c(0, 0, NA))
  }
  set.seed(97)
  before <- .Random.seed
  invisible(mc_loss(mc_premium(data.frame(pv_benefit = 2, pv_annuity = 1))))
  expect_identical(.Random.seed, before)
})

test_that("annualized input retains its coefficient after installment underflow", {
  d <- data.frame(pv_benefit = 1, pv_annuity = 1e100)
  out <- mc_loss(d, P = 1e-300, premium_unit = "annualized", k = 1e100,
                 annuity_payment = 1e-300)
  expect_equal(out$premium_per_payment, 0)
  expect_equal(out$premium_annuity_scale / 1e-100, 1, tolerance = 1e-12)
  expect_equal(out$pv_premiums, 1, tolerance = 1e-12)
  expect_lt(abs(out$L), 1e-12)
})

test_that("known zero or signed annuity values survive a lost ratio", {
  d <- data.frame(pv_benefit = c(0, 0, NA), pv_annuity = c(0, -1e-300, NA))
  out <- mc_loss(d, P = 1e300, premium_unit = "per_payment", annuity_payment = 1e-300)
  expect_identical(out$pv_premiums[1], 0)
  expect_equal(out$pv_premiums[2] / -1e300, 1, tolerance = 1e-12)
  expect_true(is.na(out$pv_premiums[3]))
  expect_equal(out$L[2] / 1e300, 1, tolerance = 1e-12)
})

test_that("partially rounded subnormal intermediates do not bias final amounts", {
  tiny <- .Machine$double.xmin * .Machine$double.eps
  p <- mc_premium(data.frame(pv_benefit = 1.5, pv_annuity = 1),
                   k = 1e300, annuity_payment = tiny)
  expect_equal(p$premium_annualized / (tiny * 1e300 * 1.5), 1, tolerance = 1e-12)
  loss <- mc_loss(data.frame(pv_benefit = 2, pv_annuity = 1),
                   P = 3 * tiny, premium_unit = "annualized", k = 2, annuity_payment = tiny)
  expect_equal(loss$premium_annuity_scale, 1.5, tolerance = 1e-12)
  expect_equal(loss$pv_premiums, 1.5, tolerance = 1e-12)
  expect_equal(loss$L, .5, tolerance = 1e-12)
})
