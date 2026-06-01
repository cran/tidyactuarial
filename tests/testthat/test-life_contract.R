test_that("life_contract creates a single-life contract with compact actuarial notation", {
  lt <- data.frame(
    x = 40:90,
    lx = round(100000 * exp(-0.018 * (0:50)^1.35))
  )
  lt$lx[nrow(lt)] <- 0

  contract <- life_contract(
    lt = lt,
    lives = "single",
    x = 40,
    i = 0.05
  )

  expect_s3_class(contract, "tidyact_life_contract")
  expect_identical(contract$lives, "single")
  expect_identical(contract$x, 40L)
  expect_null(contract$y)
  expect_identical(contract$i, 0.05)
  expect_identical(contract$i_type, "effective")
  expect_identical(contract$m, 1L)

  # Temporary internal compatibility fields during the 0.1.4 API migration.
  expect_identical(contract$mortality_table, lt)
  expect_identical(contract$age, 40L)
  expect_null(contract$age_x)
  expect_null(contract$age_y)
  expect_identical(contract$rate, 0.05)
  expect_identical(contract$rate_type, "effective")
})

test_that("life_contract creates a two-life contract with one life table", {
  lt <- data.frame(
    x = 50:110,
    lx = seq(100000, 0, length.out = 61)
  )

  contract <- life_contract(
    lt = lt,
    lives = "joint",
    x = 60,
    y = 58,
    i = 0.04,
    i_type = "effective",
    m = 1
  )

  expect_s3_class(contract, "tidyact_life_contract")
  expect_identical(contract$lives, "joint")
  expect_identical(contract$x, 60L)
  expect_identical(contract$y, 58L)
  expect_identical(contract$i, 0.04)
  expect_identical(contract$i_type, "effective")
  expect_identical(contract$m, 1L)

  # Temporary internal compatibility fields during the 0.1.4 API migration.
  expect_identical(contract$mortality_table, lt)
  expect_null(contract$age)
  expect_identical(contract$age_x, 60L)
  expect_identical(contract$age_y, 58L)
  expect_identical(contract$rate, 0.04)
  expect_identical(contract$rate_type, "effective")
})

test_that("life_contract creates a two-life contract with two independent life tables", {
  lt_m <- data.frame(
    x = 60:110,
    lx = seq(100000, 0, length.out = 51)
  )

  lt_f <- data.frame(
    x = 58:108,
    lx = seq(100000, 0, length.out = 51)
  )

  contract <- life_contract(
    lt = list(lt_m, lt_f),
    lives = "last_survivor",
    x = 60,
    y = 58,
    i = 0.05,
    i_type = "nominal_interest",
    m = 12
  )

  expect_s3_class(contract, "tidyact_life_contract")
  expect_identical(contract$lives, "last_survivor")
  expect_identical(contract$x, 60L)
  expect_identical(contract$y, 58L)
  expect_identical(contract$i, 0.05)
  expect_identical(contract$i_type, "nominal_interest")
  expect_identical(contract$m, 12L)

  expect_true(is.list(contract$lt))
  expect_length(contract$lt, 2)
  expect_identical(contract$lt[[1]], lt_m)
  expect_identical(contract$lt[[2]], lt_f)

  # Temporary internal compatibility fields during the 0.1.4 API migration.
  expect_true(is.list(contract$mortality_table))
  expect_identical(contract$age_x, 60L)
  expect_identical(contract$age_y, 58L)
})

test_that("life_contract rejects old parameter names explicitly", {
  lt <- data.frame(
    x = 40:90,
    lx = seq(100000, 0, length.out = 51)
  )

  expect_error(
    life_contract(
      mortality_table = lt,
      lives = "single",
      age = 40,
      rate = 0.05
    ),
    "Deprecated argument name"
  )

  expect_error(
    life_contract(
      lt = lt,
      lives = "single",
      x = 40,
      i = 0.05,
      rate_type = "effective"
    ),
    "Deprecated argument name"
  )
})

test_that("life_contract validates required arguments under actuarial notation", {
  lt <- data.frame(
    x = 40:90,
    lx = seq(100000, 0, length.out = 51)
  )

  expect_error(
    life_contract(
      lt = lt,
      lives = "single",
      i = 0.05
    ),
    "`x` must be provided"
  )

  expect_error(
    life_contract(
      lt = lt,
      lives = "joint",
      x = 60,
      i = 0.05
    ),
    "`x` and `y` must be provided"
  )

  expect_error(
    life_contract(
      lt = lt,
      lives = "single",
      x = 40
    ),
    "`i` must be provided"
  )

  expect_error(
    life_contract(
      lt = lt,
      lives = "single",
      x = 40,
      i = 0.05,
      i_type = "bad_type"
    ),
    "`i_type` must be one of"
  )

  expect_error(
    life_contract(
      lt = lt,
      lives = "single",
      x = 40,
      i = 0.05,
      m = 0
    ),
    "`m` must be a single positive integer"
  )
})

test_that("print.tidyact_life_contract displays compact actuarial fields", {
  lt <- data.frame(
    x = 40:90,
    lx = seq(100000, 0, length.out = 51)
  )

  contract <- life_contract(
    lt = lt,
    lives = "single",
    x = 40,
    i = 0.05
  )

  printed <- capture.output(print(contract))

  expect_true(any(grepl("<tidyact_life_contract>", printed, fixed = TRUE)))
  expect_true(any(grepl("lives:", printed, fixed = TRUE)))
  expect_true(any(grepl("x:", printed, fixed = TRUE)))
  expect_true(any(grepl("i:", printed, fixed = TRUE)))
  expect_true(any(grepl("i_type:", printed, fixed = TRUE)))
  expect_true(any(grepl("m:", printed, fixed = TRUE)))
})
