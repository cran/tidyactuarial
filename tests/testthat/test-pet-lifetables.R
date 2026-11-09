test_that("pet tables retain source-specific radix and terminal intervals", {
  z <- pet_lifetables
  expect_s3_class(z, "tbl_df")
  expect_equal(dim(z), c(57, 17))
  expect_equal(anyDuplicated(z[c("table_id", "x")]), 0)
  expect_true(all(z$qx >= 0 & z$qx <= 1))
  first <- z[z$x == 0, ]
  expect_equal(first$lx, c(100000, 100000, 30563))
  expect_equal(first$ex, c(12.69, 11.18, 11.23))
  last <- z[z$open_interval, ]
  expect_equal(last$x, c(17L, 17L, 20L))
  expect_equal(last$qx, c(.4166, .4015, 1))
  expect_equal(last$lx, c(13204, 13168, 23))
  expect_true(all(is.na(z$mx[z$country == "United Kingdom"])))
  expect_true(all(is.na(z$ax[z$country == "United States"])))
})

test_that("annual survivals agree at source precision without closing the tail", {
  groups <- split(pet_lifetables, pet_lifetables$table_id)
  expect_true(all(vapply(groups, function(g) {
    identical(g$x, 0:(nrow(g) - 1L)) && sum(g$open_interval) == 1L &&
      max(abs(1 - g$lx[-1] / head(g$lx, -1) - head(g$qx, -1))) < .00051
  }, TRUE)))
  expect_lt(max(abs(pet_lifetables$Tx / pet_lifetables$lx - pet_lifetables$ex)), .005)
  dogs <- subset(pet_lifetables, table_id == "montoya2023_us_dog", select = c(x, lx))
  expect_equal(t_px(dogs, 0, 5), .93223)
})
