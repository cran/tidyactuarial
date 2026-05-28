test_that("soa08lt loads and has valid life-table structure", {
  data(soa08lt)

  expect_s3_class(soa08lt, "data.frame")
  expect_true(all(c("x", "lx", "dx", "qx", "px") %in% names(soa08lt)))
  expect_gt(nrow(soa08lt), 0)
  expect_equal(anyDuplicated(soa08lt$x), 0L)
  expect_true(all(is.finite(soa08lt$x)))
  expect_true(all(is.finite(soa08lt$lx)))
  expect_true(all(soa08lt$lx >= 0))
  expect_true(all(soa08lt$qx >= 0 & soa08lt$qx <= 1))
  expect_true(all(soa08lt$px >= 0 & soa08lt$px <= 1))
})

test_that("soa08lt supports basic life-contingency calculations", {
  data(soa08lt)

  ax <- annuity_x(
    mortality_table = soa08lt,
    age = 65,
    rate = 0.06,
    timing = "due"
  )

  Ax <- insurance_x(
    mortality_table = soa08lt,
    age = 65,
    rate = 0.06,
    insurance_type = "whole"
  )

  expect_type(ax, "double")
  expect_type(Ax, "double")
  expect_length(ax, 1)
  expect_length(Ax, 1)
  expect_true(is.finite(ax))
  expect_true(is.finite(Ax))
  expect_gte(ax, 0)
  expect_gte(Ax, 0)
})

test_that("soa08lt satisfies the basic two-life annuity identity", {
  data(soa08lt)

  ax <- annuity_x(
    mortality_table = soa08lt,
    age = 60,
    rate = 0.06,
    timing = "due"
  )

  ay <- annuity_x(
    mortality_table = soa08lt,
    age = 65,
    rate = 0.06,
    timing = "due"
  )

  axy <- annuity_xy(
    mortality_table = soa08lt,
    age_x = 60,
    age_y = 65,
    rate = 0.06,
    cohort = "first",
    timing = "due"
  )

  alast <- annuity_xy(
    mortality_table = soa08lt,
    age_x = 60,
    age_y = 65,
    rate = 0.06,
    cohort = "last",
    timing = "due"
  )

  expect_equal(alast, ax + ay - axy, tolerance = 1e-8)
})
