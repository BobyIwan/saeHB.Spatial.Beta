test_that("moran_test works correctly with valid data", {
  # Convert weight matrix to listw object using the dataset we created
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Test default method (randomization)
  res <- moran_test(x = databeta$y, listw = W_listw)

  expect_s3_class(res, "htest")
  expect_equal(res$method, "Moran's I test under randomization (analytical)")
  expect_true(is.numeric(res$statistic))
})

test_that("moran_test handles missing values (NA) correctly", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Test with na.rm = TRUE (should output an info message and succeed)
  expect_message(
    res_na <- moran_test(x = databeta_na$y, listw = W_listw, na.rm = TRUE),
    "missing values detected and removed"
  )
  expect_s3_class(res_na, "htest")

  # Test with na.rm = FALSE (should throw an error)
  expect_error(
    moran_test(x = databeta_na$y, listw = W_listw, na.rm = FALSE),
    "Missing values \\(NA\\) detected"
  )
})

test_that("moran_test works with Monte Carlo permutation and alternative hypothesis", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Test with mc = TRUE
  res_mc <- moran_test(x = databeta$y, listw = W_listw, mc = TRUE, nsim = 99)

  expect_s3_class(res_mc, "htest")
  expect_equal(res_mc$method, "Moran's I test under Monte Carlo permutation")

  # Test non-default alternative hypothesis
  res_less <- moran_test(x = databeta$y, listw = W_listw, alternative = "less")
  expect_match(res_less$alternative, "less")

  # Test alternative = "two.sided"
  res_two <- moran_test(x = databeta$y, listw = W_listw, alternative = "two.sided")
  expect_match(res_two$alternative, "two.sided")

})

test_that("moran_test stops with invalid inputs", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Error if x is not numeric
  expect_error(
    moran_test(x = as.character(databeta$y), listw = W_listw),
    "must be a numeric vector"
  )

  # Error if listw is not of class listw (passing a regular matrix)
  expect_error(
    moran_test(x = databeta$y, listw = weight_mat),
    "must be an object of class 'listw'"
  )

  # Error if length of x does not match listw
  expect_error(
    moran_test(x = databeta$y[1:10], listw = W_listw),
    "does not match the number of spatial units"
  )

  # Error if valid data points are less than 3
  x_short <- c(1, 2, rep(NA, 34))
  expect_error(
    moran_test(x = x_short, listw = W_listw, na.rm = TRUE),
    "Not enough valid data points"
  )
})
