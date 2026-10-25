test_that("two-asset duration immunization matches liability PV and duration", {
  i <- 0.05
  L <- 1000
  t <- 5
  pv_L <- L / (1 + i)^t

  out <- immunize_duration(
    L = L,
    t = t,
    P = c(1, 1),
    D = c(3, 7),
    i = i
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$PV_L[1], pv_L, tolerance = 1e-10)
  expect_equal(out$D_L[1], 5, tolerance = 1e-12)
  expect_equal(out$PV_A[1], pv_L, tolerance = 1e-10)
  expect_equal(out$D_A[1], 5, tolerance = 1e-10)
  expect_equal(out$w, rep(pv_L / 2, 2), tolerance = 1e-10)
})

test_that("duration immunization works with non-effective input rates", {
  j12 <- 0.06
  i_eff <- (1 + j12 / 12)^12 - 1

  out <- immunize_duration(
    L = 1000,
    t = 5,
    P = c(1, 1),
    D = c(3, 7),
    i = j12,
    i_type = "nominal_interest",
    m = 12
  )

  expected_pv <- 1000 / (1 + i_eff)^5
  expect_equal(out$PV_L[1], expected_pv, tolerance = 1e-9)
  expect_equal(out$D_A[1], out$D_L[1], tolerance = 1e-10)
})

test_that("duration immunization rejects equal-duration two-asset systems", {
  expect_error(
    immunize_duration(
      L = 1000,
      t = 5,
      P = c(100, 100),
      D = c(4, 4),
      i = 0.05
    ),
    "equal duration"
  )
})

test_that("duration-convexity immunization matches all three moments", {
  i <- 0.05
  L <- 1000
  t <- 5
  pv_L <- L / (1 + i)^t

  D_assets <- c(2, 5, 8)
  C_assets <- D_assets * (D_assets + 1) / (1 + i)^2

  out <- immunize_duration_convexity(
    L = L,
    t = t,
    P = c(1, 1, 1),
    D = D_assets,
    C = C_assets,
    i = i
  )

  expected_C <- t * (t + 1) / (1 + i)^2

  expect_s3_class(out, "tbl_df")
  expect_equal(out$PV_L[1], pv_L, tolerance = 1e-10)
  expect_equal(out$D_L[1], t, tolerance = 1e-12)
  expect_equal(out$C_L[1], expected_C, tolerance = 1e-10)
  expect_equal(out$PV_A[1], out$PV_L[1], tolerance = 1e-9)
  expect_equal(out$D_A[1], out$D_L[1], tolerance = 1e-9)
  expect_equal(out$C_A[1], out$C_L[1], tolerance = 1e-9)
  expect_equal(out$w, c(0, pv_L, 0), tolerance = 1e-8)
})

test_that("duration-convexity immunization requires at least three assets", {
  expect_error(
    immunize_duration_convexity(
      L = 1000,
      t = 5,
      P = c(100, 100),
      D = c(3, 7),
      C = c(12, 56),
      i = 0.05
    ),
    "At least three assets"
  )
})

test_that("plot_immunization_gap is zero for an exactly matching asset and liability", {
  p <- plot_immunization_gap(
    L = 1000,
    t = 5,
    asset_cashflows = list(
      list(cf = 1000, t = 5)
    ),
    w = 1,
    i = 0.05,
    delta = 0.01,
    n_grid = 21
  )

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 21L)
  expect_true(all(abs(p$data$gap) < 1e-10))
})
