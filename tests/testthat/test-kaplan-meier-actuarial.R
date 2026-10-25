km_example_data <- function() {
  list(
    time = c(1, 2, 3, 4, 5, 6),
    status = c(1, 0, 1, 1, 0, 1)
  )
}

test_that("km_lifetable reproduces the Kaplan-Meier product manually", {
  dat <- km_example_data()

  out <- km_lifetable(
    time = dat$time,
    status = dat$status,
    breaks = 0:6,
    radix = 1000,
    assumption = "UDD"
  )

  km <- out$km

  expect_equal(km$time, 1:6)
  expect_equal(km$n_risk, 6:1)
  expect_equal(km$d, c(1, 0, 1, 1, 0, 1))
  expect_equal(km$censored, c(0, 1, 0, 0, 1, 0))

  expected_S <- c(
    5 / 6,
    5 / 6,
    (5 / 6) * (3 / 4),
    (5 / 6) * (3 / 4) * (2 / 3),
    (5 / 6) * (3 / 4) * (2 / 3),
    0
  )

  expect_equal(km$S, expected_S, tolerance = 1e-12)

  # Censoring alone does not create a survival jump.
  expect_equal(km$S[2], km$S[1], tolerance = 1e-12)
  expect_equal(km$S[5], km$S[4], tolerance = 1e-12)
})

test_that("km_lifetable Greenwood variance follows the actuarial formula", {
  dat <- km_example_data()

  out <- km_lifetable(
    time = dat$time,
    status = dat$status,
    breaks = 0:6
  )

  km <- out$km

  n <- c(6, 5, 4, 3, 2, 1)
  d <- c(1, 0, 1, 1, 0, 1)

  g_term <- ifelse(n > d & n > 0, d / (n * (n - d)), 0)
  expected_var <- km$S^2 * cumsum(g_term)

  expect_equal(km$varS, expected_var, tolerance = 1e-12)
  expect_equal(km$seS, sqrt(expected_var), tolerance = 1e-12)
})

test_that("km_lifetable maps KM survival into an empirical UDD life table", {
  dat <- km_example_data()

  out <- km_lifetable(
    time = dat$time,
    status = dat$status,
    breaks = 0:6,
    radix = 1000,
    assumption = "UDD"
  )

  lt <- out$lifetable

  expected_lx <- 1000 * c(
    1,
    5 / 6,
    5 / 6,
    5 / 8,
    5 / 12,
    5 / 12
  )

  expected_lx_next <- 1000 * c(
    5 / 6,
    5 / 6,
    5 / 8,
    5 / 12,
    5 / 12,
    0
  )

  expected_dx <- expected_lx - expected_lx_next
  expected_qx <- expected_dx / expected_lx
  expected_Lx <- 0.5 * (expected_lx + expected_lx_next)

  expect_s3_class(lt, "lifetable")
  expect_equal(lt$x, 0:5)
  expect_equal(lt$lx, expected_lx, tolerance = 1e-10)
  expect_equal(lt$dx, expected_dx, tolerance = 1e-10)
  expect_equal(lt$qx, expected_qx, tolerance = 1e-10)
  expect_equal(lt$px, 1 - expected_qx, tolerance = 1e-10)
  expect_equal(lt$Lx, expected_Lx, tolerance = 1e-10)

  expected_Tx <- rev(cumsum(rev(expected_Lx)))
  expect_equal(lt$Tx, expected_Tx, tolerance = 1e-10)
  expect_equal(lt$ex, expected_Tx / expected_lx, tolerance = 1e-10)

  expect_equal(attr(lt, "source"), "Kaplan-Meier")
  expect_equal(attr(lt, "frac"), "UDD")
  expect_false(attr(lt, "closed"))
})

test_that("km_lifetable handles delayed entry in the risk set", {
  out <- km_lifetable(
    time = c(1, 2, 3),
    status = c(0, 1, 1),
    entry = c(0, 1.5, 0),
    breaks = 0:3
  )

  # At t=1 the second individual has not entered yet.
  expect_equal(out$km$n_risk, c(2, 2, 1))
  expect_equal(out$km$S, c(1, 0.5, 0), tolerance = 1e-12)
})

test_that("km_lifetable uses the documented UDD, CF and Balducci exposures", {
  time <- c(0.5, 2)
  status <- c(1, 0)
  breaks <- c(0, 1)
  radix <- 1000

  udd <- km_lifetable(
    time = time,
    status = status,
    breaks = breaks,
    radix = radix,
    assumption = "UDD"
  )$lifetable

  cf <- km_lifetable(
    time = time,
    status = status,
    breaks = breaks,
    radix = radix,
    assumption = "CF"
  )$lifetable

  bal <- km_lifetable(
    time = time,
    status = status,
    breaks = breaks,
    radix = radix,
    assumption = "Balducci"
  )$lifetable

  l0 <- 1000
  l1 <- 500

  expected_udd <- 0.5 * (l0 + l1)
  expected_cf <- (l0 - l1) / log(l0 / l1)
  expected_bal <- (l0 * l1) / (l0 - l1) * log(l0 / l1)

  expect_equal(udd$Lx, expected_udd, tolerance = 1e-10)
  expect_equal(cf$Lx, expected_cf, tolerance = 1e-10)
  expect_equal(bal$Lx, expected_bal, tolerance = 1e-10)
})

test_that("km_lifetable validates censoring and entry inputs", {
  expect_error(
    km_lifetable(
      time = c(1, 2),
      status = c(1, 2),
      breaks = 0:2
    ),
    "0/1"
  )

  expect_error(
    km_lifetable(
      time = c(1, 2),
      status = c(1, 0),
      entry = c(0, 3),
      breaks = 0:2
    ),
    "entry <= time"
  )

  expect_error(
    km_lifetable(
      time = c(1, 2),
      status = c(1, 0),
      breaks = 0
    ),
    "at least two"
  )
})
