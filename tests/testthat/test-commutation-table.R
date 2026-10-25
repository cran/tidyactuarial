test_that("commutation_table reproduces classical commutation functions", {
  lt <- lifetable(
    x = 60:63,
    lx = c(100000, 90000, 72000, 36000),
    close = TRUE
  )

  i <- 0.05
  v <- 1 / (1 + i)
  x <- lt$x
  lx <- lt$lx
  dx <- lx - c(lx[-1], 0)

  Dx <- v^x * lx
  Cx <- v^(x + 1) * dx
  rev_cumsum <- function(z) rev(cumsum(rev(z)))
  Nx <- rev_cumsum(Dx)
  Mx <- rev_cumsum(Cx)
  Sx <- rev_cumsum(Nx)
  Rx <- rev_cumsum(Mx)

  out <- commutation_table(lt = lt, i = i)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("x", "lx", "dx", "v", "Dx", "Nx", "Sx", "Cx", "Mx", "Rx"))
  expect_equal(out$dx, dx, tolerance = 1e-12)
  expect_equal(out$v, rep(v, length(x)), tolerance = 1e-12)
  expect_equal(out$Dx, Dx, tolerance = 1e-12)
  expect_equal(out$Cx, Cx, tolerance = 1e-12)
  expect_equal(out$Nx, Nx, tolerance = 1e-12)
  expect_equal(out$Mx, Mx, tolerance = 1e-12)
  expect_equal(out$Sx, Sx, tolerance = 1e-12)
  expect_equal(out$Rx, Rx, tolerance = 1e-12)
})

test_that("at zero interest Mx equals lx for a closed life table", {
  lt <- lifetable(
    x = 60:63,
    lx = c(100000, 90000, 72000, 36000),
    close = TRUE
  )

  out <- commutation_table(lt = lt, i = 0)

  expect_equal(out$v, rep(1, nrow(lt)))
  expect_equal(out$Dx, lt$lx)
  expect_equal(out$Cx, lt$dx)
  expect_equal(out$Mx, lt$lx, tolerance = 1e-12)
})

test_that("commutation_table gives the same result for equivalent interest conventions", {
  lt <- lifetable(
    x = 60:63,
    lx = c(100000, 90000, 72000, 36000),
    close = TRUE
  )

  j12 <- 0.06
  i_eff <- (1 + j12 / 12)^12 - 1

  effective <- commutation_table(lt, i = i_eff)
  nominal <- commutation_table(
    lt,
    i = j12,
    i_type = "nominal_interest",
    m = 12
  )

  expect_equal(nominal$Dx, effective$Dx, tolerance = 1e-10)
  expect_equal(nominal$Cx, effective$Cx, tolerance = 1e-10)
  expect_equal(nominal$Nx, effective$Nx, tolerance = 1e-10)
  expect_equal(nominal$Mx, effective$Mx, tolerance = 1e-10)
})

test_that("commutation_table validates life table and interest inputs", {
  bad_lt <- data.frame(age = 60:62, survivors = c(1000, 900, 700))

  expect_error(
    commutation_table(bad_lt, i = 0.05),
    "must contain columns"
  )

  lt <- data.frame(x = 60:62, lx = c(1000, 900, 700))

  expect_error(
    commutation_table(lt, i = 0.05, i_type = "unknown"),
    "must be one of"
  )
})
