test_that("stop-loss premium is exact for a discrete loss distribution", {
  severity <- severity_spec("discrete", values = c(0, 100, 200),
                            probs = c(0.3, 0.4, 0.3))
  out <- stop_loss_premium(severity, retention = c(0, 150, 250))
  expect_equal(out$premium, c(100, 15, 0))
})

test_that("exponential and gamma premiums use exact formulas", {
  exp_model <- severity_spec("exponential", rate = 0.02)
  expect_equal(stop_loss_premium(exp_model, retention = 10)$premium,
               exp(-0.2) / 0.02)
  gam_model <- severity_spec("gamma", shape = 2, rate = 0.1)
  expect_equal(stop_loss_premium(gam_model, retention = 0)$premium, 20)
})

test_that("Pareto II stop-loss recognizes infinite mean", {
  heavy <- severity_spec("pareto", shape = 0.9, scale = 100)
  expect_equal(stop_loss_premium(heavy, 200)$premium, Inf)
  finite <- severity_spec("pareto", shape = 2, scale = 100)
  expect_equal(stop_loss_premium(finite, 100)$premium, 50)
})

test_that("stop-loss rejects aggregate distributions with missing tail", {
  grid <- structure(list(grid = tibble::tibble(
    loss = c(0, 50), probability = c(0.5, 0.45)),
    tail_probability = 0.05), class = "tidyact_aggregate_distribution")
  expect_error(stop_loss_premium(grid, 10), "uncomputed aggregate tail")
})

test_that("stop-loss handles simulation tibbles and invalid retentions", {
  x <- tibble::tibble(aggregate_loss = c(0, 100, 100, 300))
  expect_equal(stop_loss_premium(x, 100)$premium, 50)
  expect_error(stop_loss_premium(x, -1), "retention")
})

test_that("covered exponential severity has exact stop-loss", {
  x <- severity_spec("exponential", rate = 0.01)
  paid <- apply_coverage(x, coverage_spec(deductible = 50, limit = 200))
  out <- stop_loss_premium(paid, retention = 0)
  # E[min((X-50)+, 200)] = P(X>50)(1-exp(-.01*200))/.01
  expected <- exp(-0.5) * (1 - exp(-2)) / 0.01
  expect_equal(out$premium, expected, tolerance = 1e-5)
  expect_identical(out$method, "closed_form")
})

test_that("log-logistic with infinite mean has infinite stop-loss", {
  x <- severity_spec("loglogistic", shape = 0.8, scale = 100)
  expect_equal(stop_loss_premium(x, 100)$premium, Inf)
})
