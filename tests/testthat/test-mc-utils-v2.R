test_that(".mc_assert_character_scalar rejects empty names", {
  expect_error(
    tidyactuarial:::.mc_assert_character_scalar("", "col"),
    "nonempty"
  )

  expect_error(
    tidyactuarial:::.mc_assert_character_scalar("   ", "col"),
    "nonempty"
  )

  expect_invisible(
    tidyactuarial:::.mc_assert_character_scalar("pv_benefit", "col")
  )
})


test_that(".mc_assert_column validates both data and column names", {
  expect_error(
    tidyactuarial:::.mc_assert_column(
      data = 1:3,
      col = "x",
      arg = "col"
    ),
    "data frame"
  )

  expect_error(
    tidyactuarial:::.mc_assert_column(
      data = data.frame(x = 1:3),
      col = "y",
      arg = "col"
    ),
    "identify a column"
  )
})


test_that(".mc_effective_rate handles deprecated aliases exactly", {
  expect_equal(
    tidyactuarial:::.mc_effective_rate(
      rate = 0.05,
      interest_type = "effective"
    ),
    0.05
  )

  expect_error(
    tidyactuarial:::.mc_effective_rate(
      i = 0.05,
      rate = 0.05
    ),
    "Provide only one of `i`"
  )

  expect_error(
    tidyactuarial:::.mc_effective_rate(
      i = 0.05,
      i_type = "effective",
      interest_type = "effective"
    ),
    "Provide only one of `i_type`"
  )
})


test_that(".mc_effective_rate requires an integer conversion frequency", {
  expect_error(
    tidyactuarial:::.mc_effective_rate(
      i = 0.06,
      i_type = "nominal_interest",
      m = 1.5
    ),
    "positive integer"
  )
})


test_that(".mc_discount_factor supports deprecated rate aliases", {
  expect_equal(
    tidyactuarial:::.mc_discount_factor(
      rate = 0.05,
      interest_type = "effective"
    ),
    1 / 1.05,
    tolerance = 1e-12
  )

  expect_equal(
    tidyactuarial:::.mc_discount_factor(
      i = 0.05
    ),
    1 / 1.05,
    tolerance = 1e-12
  )

  expect_error(
    tidyactuarial:::.mc_discount_factor(
      i = 0.05,
      rate = 0.05
    ),
    "Provide only one of `i`"
  )
})


test_that(".mc_payment_times builds stable fractional grids", {
  expect_equal(
    tidyactuarial:::.mc_payment_times(
      from = 0,
      to = 1,
      by = 1 / 12
    ),
    (0:12) / 12,
    tolerance = 1e-12
  )

  expect_equal(
    tidyactuarial:::.mc_payment_times(
      from = 0.25,
      to = 1,
      by = 0.25
    ),
    c(0.25, 0.5, 0.75, 1)
  )

  expect_length(
    tidyactuarial:::.mc_payment_times(
      from = 2,
      to = 1,
      by = 0.25
    ),
    0L
  )
})


test_that(".mc_payment_times rejects invalid steps", {
  expect_error(
    tidyactuarial:::.mc_payment_times(
      from = 0,
      to = 1,
      by = 0
    ),
    "positive"
  )

  expect_error(
    tidyactuarial:::.mc_payment_times(
      from = "0",
      to = 1,
      by = 0.25
    ),
    "numeric scalars"
  )
})


test_that(".mc_quantile_names validates probabilities", {
  expect_identical(
    tidyactuarial:::.mc_quantile_names(
      c(0.025, 0.5, 0.975)
    ),
    c("q2_5", "q50", "q97_5")
  )

  expect_error(
    tidyactuarial:::.mc_quantile_names(
      c(-0.1, 0.5)
    ),
    "\\[0, 1\\]"
  )

  expect_error(
    tidyactuarial:::.mc_quantile_names(
      numeric(0)
    ),
    "nonempty"
  )
})


test_that(".mc_quantile_names rejects display-name collisions", {
  expect_error(
    tidyactuarial:::.mc_quantile_names(
      c(0.1234561, 0.1234562)
    ),
    "duplicated"
  )
})
