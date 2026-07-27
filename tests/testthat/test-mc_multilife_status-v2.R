test_that("mc_multilife_status computes joint and last-survivor lifetimes", {
  simulated <- tibble::tibble(
    sim_id = rep(1:2, each = 2),
    life_id = rep(c("x", "y"), times = 2),
    Kx = c(2, 5, 4, 1),
    Tx = c(2.4, 5.2, 4.5, 1.3)
  )

  joint <- mc_multilife_status(
    simulated,
    status = "joint"
  )

  last <- mc_multilife_status(
    simulated,
    status = "last_survivor"
  )

  expect_equal(joint$K_status, c(2L, 1L))
  expect_equal(joint$T_status, c(2.4, 1.3))
  expect_identical(joint$status, rep("joint", 2))

  expect_equal(last$K_status, c(5L, 4L))
  expect_equal(last$T_status, c(5.2, 4.5))
  expect_identical(last$status, rep("last_survivor", 2))
})


test_that("mc_multilife_status canonicalizes transitional aliases", {
  simulated <- tibble::tibble(
    sim_id = rep(1:2, each = 2),
    life_id = rep(c("x", "y"), times = 2),
    Kx = c(2, 5, 4, 1),
    Tx = c(2.4, 5.2, 4.5, 1.3)
  )

  joint <- mc_multilife_status(
    simulated,
    status = "joint_life"
  )

  last <- mc_multilife_status(
    simulated,
    status = "last"
  )

  expect_identical(joint$status, rep("joint", 2))
  expect_identical(last$status, rep("last_survivor", 2))
})


test_that("mc_multilife_status retains a custom simulation identifier", {
  simulated <- tibble::tibble(
    scenario = rep(c("a", "b"), each = 2),
    life = rep(c("x", "y"), times = 2),
    K = c(2, 5, 4, 1),
    T = c(2.4, 5.2, 4.5, 1.3)
  )

  result <- mc_multilife_status(
    simulated,
    col_sim = "scenario",
    col_life = "life",
    col_K = "K",
    col_T = "T"
  )

  expect_equal(result$scenario, c("a", "b"))
  expect_equal(result$sim_id, result$scenario)
  expect_equal(result$sim, result$scenario)
})


test_that("mc_multilife_status rejects duplicated lives within a simulation", {
  simulated <- tibble::tibble(
    sim_id = c(1, 1, 2, 2),
    life_id = c("x", "x", "x", "y"),
    Kx = c(2, 5, 4, 1),
    Tx = c(2.4, 5.2, 4.5, 1.3)
  )

  expect_error(
    mc_multilife_status(simulated),
    "exactly once"
  )
})


test_that("mc_multilife_status requires the same life set in every simulation", {
  simulated <- tibble::tibble(
    sim_id = c(1, 1, 2, 2),
    life_id = c("x", "y", "x", "z"),
    Kx = c(2, 5, 4, 1),
    Tx = c(2.4, 5.2, 4.5, 1.3)
  )

  expect_error(
    mc_multilife_status(simulated),
    "same set of life identifiers"
  )
})


test_that("mc_multilife_status validates curtate and complete lifetimes", {
  noninteger_K <- tibble::tibble(
    sim_id = rep(1:2, each = 2),
    life_id = rep(c("x", "y"), times = 2),
    Kx = c(2.5, 5, 4, 1),
    Tx = c(2.7, 5.2, 4.5, 1.3)
  )

  inconsistent_T <- tibble::tibble(
    sim_id = rep(1:2, each = 2),
    life_id = rep(c("x", "y"), times = 2),
    Kx = c(2, 5, 4, 1),
    Tx = c(3.2, 5.2, 4.5, 1.3)
  )

  expect_error(
    mc_multilife_status(noninteger_K),
    "whole-number curtate"
  )

  expect_error(
    mc_multilife_status(inconsistent_T),
    "K <= T < K \\+ 1"
  )
})


test_that("mc_multilife_status does not silently ignore a missing life column", {
  simulated <- tibble::tibble(
    sim_id = rep(1:2, each = 2),
    Kx = c(2, 5, 4, 1),
    Tx = c(2.4, 5.2, 4.5, 1.3)
  )

  expect_error(
    mc_multilife_status(simulated),
    "Missing required column"
  )

  result <- mc_multilife_status(
    simulated,
    col_life = NULL
  )

  expect_equal(result$n_lives, c(2L, 2L))
})


test_that("mc_multilife_status requires equal row counts without life identifiers", {
  simulated <- tibble::tibble(
    sim_id = c(1, 1, 2, 2, 2),
    Kx = c(2, 5, 4, 1, 3),
    Tx = c(2.4, 5.2, 4.5, 1.3, 3.7)
  )

  expect_error(
    mc_multilife_status(
      simulated,
      col_life = NULL
    ),
    "same number of life rows"
  )
})


test_that("mc_multilife_status requires at least two lives per simulation", {
  simulated <- tibble::tibble(
    sim_id = 1:2,
    life_id = c("x", "x"),
    Kx = c(2, 4),
    Tx = c(2.4, 4.5)
  )

  expect_error(
    mc_multilife_status(simulated),
    "at least two lives"
  )
})


test_that("mc_multilife_status rejects an explicit status vector", {
  simulated <- tibble::tibble(
    sim_id = rep(1:2, each = 2),
    life_id = rep(c("x", "y"), times = 2),
    Kx = c(2, 5, 4, 1),
    Tx = c(2.4, 5.2, 4.5, 1.3)
  )

  expect_error(
    mc_multilife_status(
      simulated,
      status = c("joint", "last_survivor")
    )
  )
})
