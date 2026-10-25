test_that("e_x curtate expectation equals the survival-probability sum", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  expected <- sum(c(90000, 72000, 50400, 30240) / 100000)

  expect_equal(
    e_x(lt, x = 60, type = "curtate"),
    expected,
    tolerance = 1e-12
  )
})

test_that("temporary curtate life expectancy uses only the requested horizon", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  expected <- 90000 / 100000 + 72000 / 100000

  expect_equal(
    e_x(lt, x = 60, t = 2, type = "curtate"),
    expected,
    tolerance = 1e-12
  )
})

test_that("complete temporary life expectancy under UDD equals integrated survival", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE,
    frac = "UDD"
  )

  # Integral from 0 to 2 under UDD:
  # [(l60+l61) + (l61+l62)] / (2*l60)
  expected <- ((100000 + 90000) + (90000 + 72000)) / (2 * 100000)

  expect_equal(
    e_x(lt, x = 60, t = 2, type = "complete"),
    expected,
    tolerance = 1e-12
  )
})

test_that("complete fractional life expectancy under constant force uses the closed form", {
  lt <- lifetable(
    x = 60:63,
    lx = c(100000, 80000, 60000, 30000),
    close = TRUE,
    frac = "CF"
  )

  p60 <- 0.8
  s <- 0.5
  expected <- (1 - p60^s) / (-log(p60))

  expect_equal(
    e_x(lt, x = 60, t = s, type = "complete"),
    expected,
    tolerance = 1e-12
  )
})

test_that("e_x tidy output reports horizon, type and fractional assumption", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE,
    frac = "Balducci"
  )

  out <- e_x(
    lt,
    x = c(60, 61),
    t = 2,
    type = "curtate",
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("x", "t", "type", "frac", "ex"))
  expect_equal(out$x, c(60L, 61L))
  expect_equal(out$t, c(2, 2))
  expect_true(all(out$type == "curtate"))
  expect_true(all(out$frac == "Balducci"))
})

test_that("t_Ex equals discounted survival probability", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  i <- 0.05
  t <- 2
  expected <- (1 + i)^(-t) * (72000 / 100000)

  expect_equal(
    t_Ex(lt, x = 60, t = t, i = i),
    expected,
    tolerance = 1e-12
  )
})

test_that("t_Ex standardizes nominal interest before discounting", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  j12 <- 0.06
  i_eff <- (1 + j12 / 12)^12 - 1
  expected <- (1 + i_eff)^(-2) * (72000 / 100000)

  expect_equal(
    t_Ex(
      lt,
      x = 60,
      t = 2,
      i = j12,
      i_type = "nominal_interest",
      m = 12
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("t_Ex at time zero is one and tidy output is auditable", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  expect_equal(t_Ex(lt, x = 60, t = 0, i = 0.05), 1)

  out <- t_Ex(
    lt,
    x = c(60, 61),
    t = c(1, 2),
    i = 0.05,
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_named(
    out,
    c("x", "t", "i", "i_type", "m", "i_effective", "frac", "nEx")
  )
  expect_equal(out$i_effective, rep(0.05, 2), tolerance = 1e-12)
})
