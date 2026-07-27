test_that("mc_loss uses annualized premiums with normalized k-thly annuities", {
  simulated <- tibble::tibble(
    pv_benefit = c(10000, 9000),
    pv_annuity = c(8, 7),
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12,
    premium_annualized = 1200
  )

  result <- mc_loss(simulated)

  expect_equal(result$premium_per_payment, c(100, 100))
  expect_equal(result$pv_premiums, c(9600, 8400))
  expect_equal(result$L, c(400, 600))
  expect_equal(result$loss, result$L)
})


test_that("mc_loss gives the same loss under equivalent annuity bases", {
  normalized <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 8,
    payment_per_payment = 1 / 12,
    payment_annualized = 1,
    payments_per_year = 12,
    premium_annualized = 1200
  )

  unit_per_payment <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 96,
    payment_per_payment = 1,
    payment_annualized = 12,
    payments_per_year = 12,
    premium_annualized = 1200
  )

  loss_normalized <- mc_loss(normalized)
  loss_unit <- mc_loss(unit_per_payment)

  expect_equal(loss_normalized$pv_premiums, 9600)
  expect_equal(loss_unit$pv_premiums, 9600)
  expect_equal(loss_normalized$L, loss_unit$L)
})


test_that("mc_loss accepts premiums expressed per payment", {
  simulated <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 96,
    payment_per_payment = 1,
    payments_per_year = 12,
    premium_per_payment = 100
  )

  result <- mc_loss(simulated)

  expect_equal(result$premium_annualized, 1200)
  expect_equal(result$premium_per_payment, 100)
  expect_equal(result$pv_premiums, 9600)
  expect_equal(result$L, 400)
})


test_that("mc_loss preserves the historical annuity-scale coefficient", {
  simulated <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 48,
    payment_per_payment = 0.5,
    payments_per_year = 12,
    P = 200
  )

  result <- mc_loss(simulated)

  expect_equal(result$premium_annuity_scale, 200)
  expect_equal(result$premium_per_payment, 100)
  expect_equal(result$premium_annualized, 1200)
  expect_equal(result$pv_premiums, 9600)
  expect_equal(result$L, 400)
})


test_that("mc_loss direct P can be declared annualized", {
  simulated <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 8,
    payment_per_payment = 1 / 12,
    payments_per_year = 12
  )

  result <- mc_loss(
    simulated,
    P = 1200,
    premium_unit = "annualized"
  )

  expect_equal(result$premium_per_payment, 100)
  expect_equal(result$pv_premiums, 9600)
  expect_equal(result$L, 400)
})


test_that("mc_loss requires annuity payment metadata for fractional custom inputs", {
  simulated <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 8,
    payments_per_year = 12,
    premium_annualized = 1200
  )

  expect_error(
    mc_loss(simulated),
    "cannot be inferred"
  )

  result <- mc_loss(
    simulated,
    annuity_payment = 1 / 12
  )

  expect_equal(result$pv_premiums, 9600)
})


test_that("mc_loss detects inconsistent frequency and payment metadata", {
  inconsistent_k <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 8,
    payment_per_payment = 1 / 12,
    payments_per_year = 12,
    k = 4,
    premium_annualized = 1200
  )

  expect_error(
    mc_loss(inconsistent_k),
    "inconsistent"
  )

  inconsistent_payment <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 8,
    payment_per_payment = 1 / 12,
    payment_annualized = 2,
    payments_per_year = 12,
    premium_annualized = 1200
  )

  expect_error(
    mc_loss(inconsistent_payment),
    "payment_annualized"
  )
})


test_that("mc_loss preserves missing simulated values", {
  simulated <- tibble::tibble(
    pv_benefit = c(10000, NA_real_),
    pv_annuity = c(8, 7),
    payment_per_payment = 1 / 12,
    payments_per_year = 12,
    premium_annualized = 1200
  )

  result <- mc_loss(simulated)

  expect_equal(result$L[[1]], 400)
  expect_true(is.na(result$L[[2]]))
})


test_that("mc_loss deprecated aliases use exact conflict detection", {
  simulated <- tibble::tibble(
    Z = 10000,
    Y = 8,
    P = 1200
  )

  expect_error(
    mc_loss(
      simulated,
      col_Z = "Z",
      benefit_col = "Z",
      col_Y = "Y"
    ),
    "Provide only one of `col_Z`"
  )

  result <- mc_loss(
    data = simulated,
    benefit_col = "Z",
    annuity_col = "Y",
    premium_col = "P",
    loss_col = "loss"
  )

  expect_equal(result$loss, 400)
})


test_that("mc_loss rejects loss-column collisions", {
  simulated <- tibble::tibble(
    pv_benefit = 10000,
    pv_annuity = 8,
    premium_annualized = 1200
  )

  expect_error(
    mc_loss(
      simulated,
      col_L = "premium_annualized"
    ),
    "conflicts"
  )
})
