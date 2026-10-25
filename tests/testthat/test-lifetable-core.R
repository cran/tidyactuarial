test_that("lifetable built from lx follows annual actuarial identities", {
  x <- 60:64
  lx <- c(100000, 90000, 72000, 50400, 30240)

  lt <- lifetable(x = x, lx = lx, close = TRUE, ax = 0.5)

  expected_dx <- c(10000, 18000, 21600, 20160, 30240)
  expected_qx <- expected_dx / lx
  expected_qx[length(expected_qx)] <- 1
  expected_px <- 1 - expected_qx
  expected_mx <- expected_qx / (1 - 0.5 * expected_qx)

  expect_s3_class(lt, "lifetable")
  expect_equal(lt$x, x)
  expect_equal(lt$lx, lx)
  expect_equal(lt$dx, expected_dx)
  expect_equal(lt$qx, expected_qx, tolerance = 1e-12)
  expect_equal(lt$px, expected_px, tolerance = 1e-12)
  expect_equal(lt$mx, expected_mx, tolerance = 1e-12)

  expect_equal(attr(lt, "radix"), 100000)
  expect_identical(attr(lt, "omega"), 64L)
  expect_identical(attr(lt, "type"), "ultimate")
  expect_identical(attr(lt, "frac"), "UDD")
  expect_true(attr(lt, "closed"))
  expect_equal(attr(lt, "ax"), 0.5)
})

test_that("lifetable built from qx reconstructs lx recursively", {
  x <- 60:64
  qx <- c(0.10, 0.20, 0.30, 0.40, 1.00)

  lt <- lifetable(x = x, qx = qx, radix = 100000, close = TRUE)

  expected_lx <- c(
    100000,
    100000 * 0.90,
    100000 * 0.90 * 0.80,
    100000 * 0.90 * 0.80 * 0.70,
    100000 * 0.90 * 0.80 * 0.70 * 0.60
  )

  expect_equal(lt$lx, expected_lx, tolerance = 1e-10)
  expect_equal(lt$qx, qx, tolerance = 1e-12)
  expect_equal(lt$px, 1 - qx, tolerance = 1e-12)
})

test_that("lifetable built from px agrees with construction from qx", {
  x <- 60:64
  qx <- c(0.10, 0.20, 0.30, 0.40, 1.00)
  px <- 1 - qx

  from_qx <- lifetable(x = x, qx = qx, radix = 100000, close = TRUE)
  from_px <- lifetable(x = x, px = px, radix = 100000, close = TRUE)

  expect_equal(from_px$lx, from_qx$lx, tolerance = 1e-12)
  expect_equal(from_px$qx, from_qx$qx, tolerance = 1e-12)
  expect_equal(from_px$px, from_qx$px, tolerance = 1e-12)
})

test_that("lifetable converts central death rates using the ax relationship", {
  x <- 60:64
  mx <- c(0.02, 0.03, 0.04, 0.05, 0.10)
  ax <- 0.5

  lt <- lifetable(
    x = x,
    mx = mx,
    radix = 100000,
    close = TRUE,
    ax = ax
  )

  expected_qx <- mx / (1 + ax * mx)
  expected_qx[length(expected_qx)] <- 1

  expect_equal(lt$qx, expected_qx, tolerance = 1e-12)
  expect_equal(lt$mx[-length(x)], mx[-length(x)], tolerance = 1e-12)
  expect_equal(tail(lt$mx, 1), 2, tolerance = 1e-12)
})

test_that("lifetable truncates at omega without inventing ages", {
  lt <- lifetable(
    x = 60:67,
    lx = c(100000, 98000, 95000, 91000, 86000, 80000, 73000, 65000),
    omega = 65,
    close = TRUE
  )

  expect_equal(lt$x, 60:65)
  expect_equal(nrow(lt), 6)
  expect_identical(attr(lt, "omega"), 65L)
  expect_equal(tail(lt$qx, 1), 1)
  expect_equal(tail(lt$px, 1), 0)
})

test_that("lifetable validates the fundamental annual-table inputs", {
  expect_error(
    lifetable(x = c(60, 62, 63), lx = c(1000, 900, 800)),
    "consecutive"
  )

  expect_error(
    lifetable(x = 60:62, lx = c(1000, 1100, 800)),
    "nonincreasing"
  )

  expect_error(
    lifetable(x = 60:62, qx = c(0.1, 0.2, 1)),
    "radix"
  )

  expect_error(
    lifetable(x = 60:62, lx = c(1000, 900, 700), omega = 65),
    "cannot exceed"
  )
})
