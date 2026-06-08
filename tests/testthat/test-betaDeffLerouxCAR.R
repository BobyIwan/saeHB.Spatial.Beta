test_that("Unit Testing for betaDeffLerouxCAR: Execution and Output Format", {

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- betaDeffLerouxCAR(
      formula = y ~ x1 + x2,
      DEFF = "deff",
      n_i = "n_i",
      proxmat = adjacency_mat,
      data = dataBeta,
      plot = TRUE,
      keep.fit = TRUE
    )
    dev.off()

    expect_true(is.list(res_sampled))
    expect_true(all(c("Est", "refVar", "randeff", "coefficient", "fit") %in% names(res_sampled)))
    expect_equal(nrow(res_sampled$Est), nrow(dataBeta))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- betaDeffLerouxCAR(
      formula = y ~ x1 + x2,
      DEFF = "deff",
      n_i = "n_i",
      proxmat = adjacency_mat,
      data = dataBeta_NA,
      plot = FALSE
    )

    expect_true(is.list(res_nonsampled))
    expect_equal(nrow(res_nonsampled$Est), nrow(dataBeta_NA))

    # Case 3: Execution works with explicitly defined coef and var.coef
    res_coef <- betaDeffLerouxCAR(
      formula = y ~ x1 + x2,
      DEFF = "deff",
      n_i = "n_i",
      proxmat = adjacency_mat,
      data = dataBeta,
      coef = c(0, 0, 0),
      var.coef = c(1, 1, 1),
      plot = FALSE
    )
    expect_true(is.list(res_coef))
  })
})

test_that("Unit Testing for betaDeffLerouxCAR: Error Handling", {
  # Case 4: Response variable (y) not between 0 and 1
  data_invalid_y <- dataBeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    betaDeffLerouxCAR(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 5: Auxiliary variable (X) contains NA values
  data_invalid_x <- dataBeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    betaDeffLerouxCAR(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 6: Iteration update is less than 3
  expect_error(
    betaDeffLerouxCAR(y ~ x1 + x2, "deff", "n_i", adjacency_mat, dataBeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 7: Sample size (n_i) is less than or equal to DEFF
  data_invalid_n <- dataBeta
  data_invalid_n$n_i[1] <- 1.5
  data_invalid_n$deff[1] <- 2.0
  expect_error(
    betaDeffLerouxCAR(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_invalid_n, plot = FALSE),
    "There is at least one sampled area where n_i <= DEFF. Effective sample size must be > 1"
  )

  # Case 8: Proximity matrix dimension mismatch
  wrong_W <- adjacency_mat[1:10, 1:10]
  expect_error(
    betaDeffLerouxCAR(y ~ x1 + x2, "deff", "n_i", wrong_W, dataBeta, plot = FALSE),
    "Proximity matrix must be N x N"
  )

  # Case 9: Proximity matrix contains NA
  na_W <- adjacency_mat
  na_W[1, 2] <- NA
  expect_error(
    betaDeffLerouxCAR(y ~ x1 + x2, "deff", "n_i", na_W, dataBeta, plot = FALSE),
    "Proximity matrix contains NA"
  )

  # Case 10: Formula without predictor
  expect_error(
    betaDeffLerouxCAR(y ~ 1, "deff", "n_i", adjacency_mat, dataBeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
})
