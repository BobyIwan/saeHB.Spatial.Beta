test_that("Unit Testing for betaSAR: Execution and Output Format", {
  skip_on_cran()

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- betaSAR(
      formula = y ~ x1 + x2,
      proxmat = weight_mat,
      data = dataBeta,

      iter.mcmc = 100,
      burn.in = 50
    )
    dev.off()
    expect_true(is.list(res_sampled))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- betaSAR(
      formula = y ~ x1 + x2,
      proxmat = weight_mat,
      data = dataBeta_NA,

      iter.mcmc = 100,
      burn.in = 50,

      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
  })
})

test_that("Unit Testing for betaSAR: Error Handling", {
  # Case 3: Response variable (y) not between 0 and 1
  data_invalid_y <- dataBeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    betaSAR(y ~ x1 + x2, proxmat = weight_mat, data = data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 4: Auxiliary variable (X) contains NA values
  data_invalid_x <- dataBeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    betaSAR(y ~ x1 + x2, proxmat = weight_mat, data = data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 5: Iteration update is less than 3
  expect_error(
    betaSAR(y ~ x1 + x2, proxmat = weight_mat, data = dataBeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 6: Proximity matrix dimension mismatch
  wrong_W <- weight_mat[1:10, 1:10]
  expect_error(
    betaSAR(y ~ x1 + x2, proxmat = wrong_W, data = dataBeta, plot = FALSE),
    "Proximity matrix must be N x N"
  )

  # Case 7: Proximity matrix contains NA
  na_W <- weight_mat
  na_W[1, 2] <- NA
  expect_error(
    betaSAR(y ~ x1 + x2, proxmat = na_W, data = dataBeta, plot = FALSE),
    "Proximity matrix contains NA"
  )

  # Case 8: Formula without predictor
  expect_error(
    betaSAR(y ~ 1, proxmat = weight_mat, data = dataBeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
})
