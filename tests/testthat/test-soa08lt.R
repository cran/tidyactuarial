test_that("soa08lt has the expected actuarial structure", {
  data("soa08lt", envir = environment())

  expect_s3_class(soa08lt, "data.frame")

  expect_true(all(c("x", "lx", "dx", "qx", "px") %in% names(soa08lt)))

  expect_true(is.numeric(soa08lt$x))
  expect_true(is.numeric(soa08lt$lx))
  expect_true(is.numeric(soa08lt$dx))
  expect_true(is.numeric(soa08lt$qx))
  expect_true(is.numeric(soa08lt$px))

  expect_false(anyNA(soa08lt$x))
  expect_false(anyNA(soa08lt$lx))
  expect_false(anyNA(soa08lt$qx))
  expect_false(anyNA(soa08lt$px))

  expect_true(all(soa08lt$x == round(soa08lt$x)))
  expect_true(all(soa08lt$lx >= 0))
  expect_true(all(soa08lt$qx >= 0 & soa08lt$qx <= 1))
  expect_true(all(soa08lt$px >= 0 & soa08lt$px <= 1))

  expect_equal(soa08lt$px, 1 - soa08lt$qx, tolerance = 1e-10)
})


test_that("soa08lt supports basic single-life calculations", {
  data("soa08lt", envir = environment())

  ax <- annuity_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    timing = "due"
  )

  Ax <- insurance_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    type = "whole"
  )

  expect_true(is.numeric(ax))
  expect_true(is.numeric(Ax))

  expect_length(ax, 1)
  expect_length(Ax, 1)

  expect_true(is.finite(ax))
  expect_true(is.finite(Ax))

  expect_gt(ax, 0)
  expect_gte(Ax, 0)
  expect_lte(Ax, 1)
})


test_that("soa08lt supports basic term and endowment insurance calculations", {
  data("soa08lt", envir = environment())

  term_A <- insurance_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    type = "term",
    n = 10
  )

  endowment_A <- insurance_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    type = "endowment",
    n = 10
  )

  expect_true(is.numeric(term_A))
  expect_true(is.numeric(endowment_A))

  expect_length(term_A, 1)
  expect_length(endowment_A, 1)

  expect_true(is.finite(term_A))
  expect_true(is.finite(endowment_A))

  expect_gte(term_A, 0)
  expect_gte(endowment_A, term_A)
  expect_lte(endowment_A, 1)
})


test_that("soa08lt satisfies the basic two-life annuity identity", {
  data("soa08lt", envir = environment())

  ax <- annuity_x(
    lt = soa08lt,
    x = 60,
    i = 0.06,
    timing = "due"
  )

  ay <- annuity_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    timing = "due"
  )

  axy_joint <- annuity_xy(
    lt = soa08lt,
    x = 60,
    y = 65,
    i = 0.06,
    status = "joint",
    timing = "due"
  )

  axy_last <- annuity_xy(
    lt = soa08lt,
    x = 60,
    y = 65,
    i = 0.06,
    status = "last",
    timing = "due"
  )

  expect_true(is.numeric(axy_joint))
  expect_true(is.numeric(axy_last))

  expect_length(axy_joint, 1)
  expect_length(axy_last, 1)

  expect_true(is.finite(axy_joint))
  expect_true(is.finite(axy_last))

  expect_gt(axy_joint, 0)
  expect_gt(axy_last, 0)
  expect_gte(axy_last, axy_joint)

  # Classical independent-lives identity:
  # a_last = a_x + a_y - a_joint
  expect_equal(
    axy_last,
    ax + ay - axy_joint,
    tolerance = 1e-8
  )
})


test_that("soa08lt works with tidy outputs where available", {
  data("soa08lt", envir = environment())

  ax_tbl <- annuity_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    timing = "due",
    tidy = TRUE
  )

  Ax_tbl <- insurance_x(
    lt = soa08lt,
    x = 65,
    i = 0.06,
    type = "whole",
    tidy = TRUE
  )

  expect_s3_class(ax_tbl, "data.frame")
  expect_s3_class(Ax_tbl, "data.frame")

  expect_equal(nrow(ax_tbl), 1)
  expect_equal(nrow(Ax_tbl), 1)

  expect_true("apv" %in% names(ax_tbl))
  expect_true("apv" %in% names(Ax_tbl))

  expect_true(is.finite(ax_tbl$apv[[1]]))
  expect_true(is.finite(Ax_tbl$apv[[1]]))
})
