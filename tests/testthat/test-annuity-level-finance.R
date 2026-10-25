test_that("a_angle matches direct present-value sums for immediate and due annuities", {
  i <- 0.05
  n <- 10
  v <- 1 / (1 + i)

  expected_immediate <- sum(v^(seq_len(n)))
  expected_due <- sum(v^(0:(n - 1)))

  expect_equal(
    a_angle(n = n, i = i, timing = "immediate"),
    expected_immediate,
    tolerance = 1e-12
  )

  expect_equal(
    a_angle(n = n, i = i, timing = "due"),
    expected_due,
    tolerance = 1e-12
  )
})

test_that("a_angle satisfies the standard due-immediate identity", {
  i <- 0.07
  n <- 12

  immediate <- a_angle(n = n, i = i, timing = "immediate")
  due <- a_angle(n = n, i = i, timing = "due")

  expect_equal(due, (1 + i) * immediate, tolerance = 1e-12)
})

test_that("a_angle applies deferment as additional discounting", {
  i <- 0.04
  n <- 8
  h <- 3

  base <- a_angle(n = n, i = i, timing = "immediate")
  deferred <- a_angle(n = n, i = i, h = h, timing = "immediate")

  expect_equal(deferred, base * (1 + i)^(-h), tolerance = 1e-12)
})

test_that("a_angle handles zero interest and continuous force correctly", {
  expect_equal(
    a_angle(n = 6, i = 0, timing = "immediate"),
    6,
    tolerance = 1e-12
  )

  delta <- 0.035
  n <- 9
  expected_continuous <- (1 - exp(-delta * n)) / delta

  expect_equal(
    a_angle(n = n, i = delta, i_type = "force", timing = "continuous"),
    expected_continuous,
    tolerance = 1e-10
  )
})

test_that("a_angle supports standard perpetuity identities", {
  i <- 0.06

  expect_equal(
    a_angle(n = NULL, i = i, perpetuity = TRUE, timing = "immediate"),
    1 / i,
    tolerance = 1e-12
  )

  expect_equal(
    a_angle(n = NULL, i = i, perpetuity = TRUE, timing = "due"),
    (1 + i) / i,
    tolerance = 1e-12
  )
})

test_that("a_angle tidy output reports the monetary present value", {
  out <- a_angle(
    n = 10,
    i = 0.05,
    payment = 1250,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$present_value, out$payment * out$annuity_factor, tolerance = 1e-12)
})

test_that("s_angle matches direct accumulated-value sums for immediate and due annuities", {
  i <- 0.05
  n <- 10

  expected_immediate <- sum((1 + i)^(0:(n - 1)))
  expected_due <- sum((1 + i)^(1:n))

  expect_equal(
    s_angle(n = n, i = i, timing = "immediate"),
    expected_immediate,
    tolerance = 1e-12
  )

  expect_equal(
    s_angle(n = n, i = i, timing = "due"),
    expected_due,
    tolerance = 1e-12
  )
})

test_that("s_angle and a_angle satisfy the standard accumulation identity", {
  i <- 0.055
  n <- 14

  expect_equal(
    s_angle(n = n, i = i, timing = "immediate"),
    (1 + i)^n * a_angle(n = n, i = i, timing = "immediate"),
    tolerance = 1e-11
  )
})

test_that("s_angle is invariant to pure deferment under its terminal-horizon convention", {
  base <- s_angle(n = 8, i = 0.04, h = 0)
  shifted <- s_angle(n = 8, i = 0.04, h = 5)

  expect_equal(shifted, base, tolerance = 1e-12)
})

test_that("s_angle handles zero interest and continuous force correctly", {
  expect_equal(s_angle(n = 7, i = 0), 7, tolerance = 1e-12)

  delta <- 0.03
  n <- 11
  expected_continuous <- (exp(delta * n) - 1) / delta

  expect_equal(
    s_angle(n = n, i = delta, i_type = "force", timing = "continuous"),
    expected_continuous,
    tolerance = 1e-10
  )
})

test_that("s_angle tidy output reports the monetary future value", {
  out <- s_angle(
    n = 10,
    i = 0.05,
    payment = 900,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$future_value, out$payment * out$accumulation_factor, tolerance = 1e-12)
})
