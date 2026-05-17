test_that("sample datasets load correctly", {
  data(mortality_world_sample_2023)
  data(mortality_world_sample_2015_2023)
  data(mortality_colombia_tables)
  data(cash_flows_sample)
  data(bonds_sample)
  data(loans_sample)
  data(multiple_decrement_sample)

  expect_s3_class(mortality_world_sample_2023, "data.frame")
  expect_s3_class(mortality_world_sample_2015_2023, "data.frame")
  expect_s3_class(mortality_colombia_tables, "data.frame")
  expect_s3_class(cash_flows_sample, "data.frame")
  expect_s3_class(bonds_sample, "data.frame")
  expect_s3_class(loans_sample, "data.frame")
  expect_s3_class(multiple_decrement_sample, "data.frame")

  expect_gt(nrow(mortality_world_sample_2023), 0)
  expect_gt(nrow(mortality_world_sample_2015_2023), 0)
  expect_gt(nrow(mortality_colombia_tables), 0)
})

test_that("mortality_law_table runs for Gompertz", {
  tab <- mortality_law_table("Gompertz", 0, 110, B = 1e-5, c = 1.08)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("x", "qx", "px", "lx", "dx") %in% names(tab)))
  expect_true(all(tab$qx >= 0 & tab$qx <= 1))
})

test_that("premium_xy returns a finite numeric value", {
  lt <- data.frame(
    x = 60:110,
    lx = seq(100000, 0, length.out = 51)
  )

  p <- premium_xy(
    lt = lt,
    x = 60,
    y = 62,
    i = 0.05,
    product = "term",
    cohort = "first",
    n = 5,
    n_prem = 3,
    benefit = 100000
  )

  expect_type(p, "double")
  expect_length(p, 1)
  expect_true(is.finite(p))
  expect_gte(p, 0)
})
