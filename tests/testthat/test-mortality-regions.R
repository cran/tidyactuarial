test_that("WPP 2023 tables preserve the official selection and age intervals", {
  z <- mortality_regions
  expect_s3_class(z, "tbl_df")
  expect_equal(dim(z), c(28785L, 21L))
  expect_identical(unique(z$year), 2023L)
  expect_setequal(unique(z$sex), c("both", "male", "female"))
  expect_length(unique(z$location_code[z$level == "country"]), 67L)
  expect_length(unique(z$location_code[z$level != "country"]), 28L)
  expect_identical(unique(z$location_code[z$location == "Latin America and the Caribbean"]), 904L)
  expect_identical(unique(z$location_code[z$location == "Australia/New Zealand"]), 927L)
  expect_equal(anyDuplicated(z[c("location_code", "year", "sex", "x")]), 0L)
  expect_true(all(z$qx >= 0 & z$qx <= 1 & z$px >= 0 & z$px <= 1))
  expect_lt(max(abs(z$qx + z$px - 1)), 1e-8)
  open <- z[z$open_interval, ]
  expect_identical(z$open_interval, z$x == 100L)
  expect_equal(nrow(open), 285L)
  expect_true(all(open$age_interval == "100+" & open$age_span == -1L &
                    open$qx == 1 & open$px == 0 & open$lx > 0))
  expect_equal(sum(abs(open$dx - open$lx) > .002), 84L)
  expect_lt(max(abs(z$Tx / z$lx - z$ex)), .0002)
})

test_that("annual WPP identities hold at the precision of the original files", {
  groups <- split(mortality_regions, list(mortality_regions$location_code,
                                        mortality_regions$sex), drop = TRUE)
  validate <- function(g) {
    a <- g[!g$open_interval, ]
    identical(g$x, 0:100) && all(a$age_span == 1L) &&
      g$lx[1L] == 1e5 && all(diff(g$lx) <= 0) &&
      max(abs(-diff(g$lx) - a$dx)) < .0016 &&
      max(abs(-diff(g$Tx) - a$Lx)) < .0016 &&
      max(abs(a$lx * a$qx - a$dx)) < .002 &&
      max(abs(a$mx * a$Lx - a$dx)) < .002
  }
  expect_length(groups, 285L)
  expect_true(all(vapply(groups, validate, TRUE)))
})

test_that("selected WPP observations retain base-1 qx and source life expectancy", {
  z <- subset(mortality_regions, location_code == 170L & x == 65L)
  z <- z[match(c("both", "female", "male"), z$sex), ]
  expect_equal(z$qx, c(.01115004, .00880358, .01379633), tolerance = 1e-12)
  expect_equal(z$lx, c(85491.226, 89108.544, 81889.090))
  expect_equal(z$ex, c(18.3787, 19.6962, 16.8579))
  col <- subset(mortality_regions, location_code == 170L & sex == "both",
                select = c(x, lx))
  expect_equal(t_px(col, x = 65, t = 5), 79847.433 / 85491.226)
})
