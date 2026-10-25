test_that("two-life consumers reject incompatible annual mortality tables", {
  valid <- data.frame(x = 60:63, lx = c(100, 80, 40, 0))
  bad <- list(
    data.frame(x = c(60, 62, 63), lx = c(100, 40, 0)),
    data.frame(x = c(60, 60, 61), lx = c(100, 80, 40)),
    data.frame(x = c(60, 61.5, 62), lx = c(100, 80, 40)),
    data.frame(x = 60:62, lx = c(100, 110, 0)),
    data.frame(x = 60:62, lx = c(100, -1, 0)),
    data.frame(x = 60:62, lx = c(100, Inf, 0)),
    data.frame(x = numeric(), lx = numeric()),
    transform(valid, qx = c(.9, .5, 1, NA)),
    transform(valid, px = c(.9, .5, 0, NA)),
    transform(valid, qx = c(.2, .5, 1, 1.1)),
    transform(valid, px = c(.8, .5, 0, NA), qx = c(.2, .4, 1, NA))
  )
  for (lt in bad) {
    expect_error(t_pxy(lt, 60, 60, .5))
    expect_error(t_pxy(list(valid, lt), 60, 60, .5))
    expect_error(e_xy(lt, 60, 60, t = 1))
  }
})

test_that("two-life age coercion and controls fail explicitly without warnings", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 40, 0))
  for (age in c(Inf, -Inf, 3e9)) {
    expect_error(t_pxy(lt, age, 60, 1), "age|finite|integer")
    expect_error(e_xy(lt, age, 60, t = 1), "age|finite|integer")
  }
  expect_error(t_pxy(lt, 60, 60, Inf), "finite")
  for (control in list(NA, 1, logical())) {
    expect_error(e_xy(lt, 60, 60, t = 1, check = control), "check")
    expect_error(e_xy(lt, 60, 60, t = 1, tidy = control), "tidy")
  }
  for (tolerance in list(-1, Inf, NA, numeric()))
    expect_error(e_xy(lt, 60, 60, t = 1, tol = tolerance), "tol")
})

test_that("valid annual tables retain independent survival and expectation identities", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 40, 0),
                   px = c(.8, .5, 0, NA), qx = c(.2, .5, 1, NA))
  shuffled <- lt[c(3, 1, 4, 2), ]
  expect_equal(t_pxy(shuffled, 60, 61, 1), .8*.5)
  expect_equal(t_pxy(list(shuffled, lt), 60, 61, 1, status = "last"), .9)
  expect_equal(e_xy(shuffled, 60, 61, t = 3), .4)
  expect_equal(e_xy(shuffled, 60, 61, t = 3, cohort = "last"), 1.3)
  # Linear UDD marginals on the first year: (1-.2u)(1-.5u).
  expect_equal(e_xy(lt, 60, 61, t = 1, type = "complete"), 1-.7/2+.1/3)
  expect_equal(e_xy(lt, 60, 61, t = 1, type = "complete", cohort = "last"),
               1-.1/3)
  expect_true(is.na(e_xy(lt, 59, 60, t = 0)))
  expect_true(is.na(e_xy(lt, 63, 60, t = 0)))
  expect_equal(t_pxy(lt, 59, 60, 0), 1)
  expect_equal(t_pxy(lt, 60, 60, 10), 0)
  expect_equal(t_pxy(lt, 60, 60, 3), 0)
  expect_equal(e_xy(lt, 60, 61, t = 1, check = FALSE), .4)
})
