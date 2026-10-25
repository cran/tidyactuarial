test_that("t_qx equals one minus t_px for integer horizons", {
  lt <- lifetable(
    x = 60:65,
    lx = c(100000, 90000, 72000, 50400, 30240, 15120),
    close = TRUE
  )

  expected_2qx <- 1 - 72000 / 100000

  expect_equal(t_qx(lt, x = 60, t = 2), expected_2qx, tolerance = 1e-12)
  expect_equal(
    t_qx(lt, x = 60, t = 2),
    1 - t_px(lt, x = 60, t = 2),
    tolerance = 1e-12
  )
})

test_that("t_qx respects UDD, constant force and Balducci fractional assumptions", {
  lt <- lifetable(
    x = 60:63,
    lx = c(100000, 80000, 60000, 30000),
    close = TRUE
  )

  q60 <- 0.20
  p60 <- 0.80
  s <- 0.5

  expected_udd <- 1 - (1 - s * q60)
  expected_cf <- 1 - p60^s
  expected_balducci <- 1 - p60 / (1 - (1 - s) * q60)

  expect_equal(t_qx(lt, 60, s, frac = "UDD"), expected_udd, tolerance = 1e-12)
  expect_equal(t_qx(lt, 60, s, frac = "CF"), expected_cf, tolerance = 1e-12)
  expect_equal(t_qx(lt, 60, s, frac = "CML"), expected_cf, tolerance = 1e-12)
  expect_equal(
    t_qx(lt, 60, s, frac = "Balducci"),
    expected_balducci,
    tolerance = 1e-12
  )
})

test_that("t_qx combines integer and fractional survival correctly", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  # 1.5 q_60 under UDD:
  # 1.5 p_60 = 1 p_60 * 0.5 p_61
  #            = (90000/100000) * (1 - 0.5*0.20)
  expected_survival <- (90000 / 100000) * (1 - 0.5 * 0.20)
  expected_death <- 1 - expected_survival

  expect_equal(
    t_qx(lt, x = 60, t = 1.5, frac = "UDD"),
    expected_death,
    tolerance = 1e-12
  )
})

test_that("t_qx returns zero at t zero and one beyond the limiting age", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE
  )

  expect_equal(t_qx(lt, 60, 0), 0)
  expect_equal(t_qx(lt, 63, 2), 1)
})

test_that("t_qx tidy output preserves the actuarial inputs", {
  lt <- lifetable(
    x = 60:64,
    lx = c(100000, 90000, 72000, 50400, 30240),
    close = TRUE,
    frac = "CF"
  )

  out <- t_qx(
    lt,
    x = c(60, 61),
    t = c(0.5, 1),
    tidy = TRUE
  )

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("x", "t", "frac", "tqx"))
  expect_equal(out$x, c(60L, 61L))
  expect_equal(out$t, c(0.5, 1))
  expect_true(all(out$frac == "CF"))
  expect_equal(
    out$tqx,
    1 - t_px(lt, x = c(60, 61), t = c(0.5, 1)),
    tolerance = 1e-12
  )
})
