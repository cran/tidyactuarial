test_that("consolidated contract components preserve single-life behavior", {
  lt <- data.frame(x = 40:90,
                   lx = round(100000 * exp(-0.012 * (0:50)^1.35)))
  single <- life_contract(lt = lt, lives = "single", x = 50, i = 0.04) |>
    add_insurance(type = "term", benefit = 100000, n = 15) |>
    add_premium_schedule(k = 12, n_prem = 10, timing = "due")
  expect_identical(single$insurance$type, "term")
  expect_identical(single$premium_schedule$k, 12L)
  expect_null(single$insurance$status)
})

test_that("consolidated contract components preserve joint and last-survivor", {
  lt <- data.frame(x = 40:90,
                   lx = round(100000 * exp(-0.012 * (0:50)^1.35)))
  joint <- life_contract(lt = lt, lives = "joint", x = 50, y = 52, i = 0.04) |>
    add_insurance(type = "term", benefit = 50000, n = 10) |>
    add_premium_schedule(k = 4, n_prem = 5)
  last <- life_contract(lt = lt, lives = "last_survivor", x = 50, y = 52,
                        i = 0.04) |>
    add_insurance(type = "pure_endowment", benefit = 30000, n = 8)
  expect_identical(joint$insurance$status, "joint")
  expect_identical(last$insurance$status, "last")
  expect_identical(joint$premium_schedule$k, 4L)
})
