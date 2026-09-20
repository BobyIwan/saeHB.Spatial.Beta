test_that("Unit Testing for moran_test: Execution and Output Format", {
  # Convert weight matrix to listw object using the dataset
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Case 1: Randomisation approach (default, mc = FALSE) executes successfully
  res_analytical <- moran_test(x = databeta$y, listw = W_listw)
  expect_s3_class(res_analytical, "htest")
  expect_true(is.numeric(res_analytical$statistic))

  # Case 2: Monte Carlo approach (mc = TRUE) executes successfully
  res_mc <- moran_test(x = databeta$y, listw = W_listw, mc = TRUE, nsim = 99)
  expect_s3_class(res_mc, "htest")

  # Case 3: Alternative hypotheses execute correctly
  res_less <- moran_test(x = databeta$y, listw = W_listw, alternative = "less")
  expect_match(res_less$alternative, "less")

  res_two <- moran_test(x = databeta$y, listw = W_listw, alternative = "two.sided")
  expect_match(res_two$alternative, "two.sided")

  # Case 4: Data with missing values (NA) executes successfully when na.rm = TRUE
  expect_message(
    res_na <- moran_test(x = databeta_na$y, listw = W_listw, na.rm = TRUE),
    "missing values detected and removed"
  )
  expect_s3_class(res_na, "htest")
})

test_that("Unit Testing for moran_test: Error Handling", {
  W_listw <- spdep::mat2listw(weight_mat, style = "W", zero.policy = TRUE)

  # Case 5: Response variable (x) is not numeric
  expect_error(
    moran_test(x = as.character(databeta$y), listw = W_listw),
    "must be a numeric vector"
  )

  # Case 6: Spatial weights (listw) is not a valid listw object
  expect_error(
    moran_test(x = databeta$y, listw = weight_mat),
    "must be an object of class 'listw'"
  )

  # Case 7: Length mismatch between x and listw
  expect_error(
    moran_test(x = databeta$y[1:10], listw = W_listw),
    "does not match the number of spatial units"
  )

  # Case 8: Missing values detected with na.rm = FALSE
  expect_error(
    moran_test(x = databeta_na$y, listw = W_listw, na.rm = FALSE),
    "Missing values \\(NA\\) detected"
  )

  # Case 9: Not enough valid data points (< 3) after removing NA
  x_short <- c(1, 2, rep(NA, length(databeta$y) - 2))
  expect_error(
    moran_test(x = x_short, listw = W_listw, na.rm = TRUE),
    "Not enough valid data points"
  )

  # Case 10: Invalid nsim argument when mc = TRUE
  expect_error(
    moran_test(x = databeta$y, listw = W_listw, mc = TRUE, nsim = 0),
    "must be a single positive integer"
  )
  expect_error(
    moran_test(x = databeta$y, listw = W_listw, mc = TRUE, nsim = 99.5),
    "must be a single positive integer"
  )
})
