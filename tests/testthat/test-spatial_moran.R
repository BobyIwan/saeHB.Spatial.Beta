test_that("spatial_moran works correctly with valid data", {
  # Convert weight matrix to listw object using the dataset we created
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Test default method (randomization)
  res <- spatial_moran(x = dataBeta$y, listw = W_listw)

  expect_s3_class(res, "htest")
  expect_equal(res$method, "Moran's I test under randomization")
  expect_true(is.numeric(res$statistic))
})

test_that("spatial_moran handles missing values (NA) correctly", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Test with na.rm = TRUE (should output an info message and succeed)
  expect_message(
    res_na <- spatial_moran(x = dataBeta_NA$y, listw = W_listw, na.rm = TRUE),
    "missing values detected and removed"
  )
  expect_s3_class(res_na, "htest")

  # Test with na.rm = FALSE (should throw an error)
  expect_error(
    spatial_moran(x = dataBeta_NA$y, listw = W_listw, na.rm = FALSE),
    "Missing values \\(NA\\) detected"
  )
})

test_that("spatial_moran works with Monte Carlo permutation and alternative hypothesis", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Test with mc = TRUE
  res_mc <- spatial_moran(x = dataBeta$y, listw = W_listw, mc = TRUE, nsim = 99)

  expect_s3_class(res_mc, "htest")
  expect_equal(res_mc$method, "Moran's I test under Monte Carlo permutation")

  # Test non-default alternative hypothesis
  res_less <- spatial_moran(x = dataBeta$y, listw = W_listw, alternative = "less")
  expect_match(res_less$alternative, "less")
})

test_that("spatial_moran stops with invalid inputs", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Error if x is not numeric
  expect_error(
    spatial_moran(x = as.character(dataBeta$y), listw = W_listw),
    "must be a numeric vector"
  )

  # Error if listw is not of class listw (passing a regular matrix)
  expect_error(
    spatial_moran(x = dataBeta$y, listw = weight_mat),
    "must be an object of class 'listw'"
  )

  # Error if valid data points are less than 3
  x_short <- c(1, 2, rep(NA, 98))
  expect_error(
    spatial_moran(x = x_short, listw = W_listw, na.rm = TRUE),
    "Not enough valid data points"
  )
})
