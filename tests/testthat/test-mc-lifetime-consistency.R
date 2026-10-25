test_that("Monte Carlo valuations reject incompatible lifetime trajectories", {
  valuators <- list(
    function(data) mc_insurance(data, i = .05),
    function(data) mc_annuity(data, i = .05),
    function(data) mc_reserve(data, t = 0, i = .05, P = 1))
  invalid <- list(
    data.frame(Kx = 1.5, Tx = 1.7),
    data.frame(Kx = 1, Tx = 4.2),
    data.frame(Kx = 2, Tx = 1.9),
    data.frame(Kx = 1, Tx = 2),
    data.frame(Kx = Inf, Tx = NA_real_),
    data.frame(Kx = 1, Tx = Inf))
  # K=floor(T): K must be a finite nonnegative integer and, when
  # both are observed, the same life must satisfy K <= T < K+1.
  for (value in valuators) {
    for (data in invalid) expect_error(value(data), "lifetime|lifetimes")
  }
})

test_that("valid and missing lifetime rows retain their existing meaning", {
  data <- data.frame(Kx = c(1, NA_real_), Tx = c(1.2, NA_real_))
  insured <- mc_insurance(data, i = 0, benefit = 100)
  annuity <- mc_annuity(data, i = 0, payment = 1, timing = "due")
  reserve <- mc_reserve(data, t = 0, i = 0, P = 10,
    premium_unit = "annualized", benefit = 100)
  expect_equal(insured$pv_benefit, c(100, NA_real_))
  expect_equal(annuity$pv_annuity, c(2, NA_real_))
  expect_equal(reserve$L_t, c(80, NA_real_))
  # A missing curtate value does not invalidate a known complete lifetime
  # for an operation that uses only the latter.
  complete <- mc_insurance(data.frame(Kx = NA_real_, Tx = 1.2),
    i = 0, benefit = 100, timing = "moment_of_death")
  expect_equal(complete$pv_benefit, 100)
  certain <- mc_annuity(data.frame(Kx = 1.5, Tx = 99), i = 0,
                        type = "certain", n = 2)
  expect_equal(certain$pv_annuity, 2)
})
