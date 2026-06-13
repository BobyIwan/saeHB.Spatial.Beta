test_that("Unit Testing for betadeff_sar: Execution and Output Format", {
  skip_on_cran()

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- betadeff_sar(
      formula = y ~ x1 + x2,
      deff = "deff",
      n_i = "n_i",
      proxmat = weight_mat,
      data = databeta,

      iter.mcmc = 100,
      burn.in = 50,
      quiet = TRUE
    )
    dev.off()
    expect_true(is.list(res_sampled))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- betadeff_sar(
      formula = y ~ x1 + x2,
      deff = "deff",
      n_i = "n_i",
      proxmat = weight_mat,
      data = databeta_na,

      iter.mcmc = 100,
      burn.in = 50,
      n.adapt = 500,

      quiet = TRUE,
      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
  })
})

test_that("Unit Testing for betadeff_sar: Error Handling", {
  # Case 3: Response variable (y) not between 0 and 1
  data_invalid_y <- databeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    betadeff_sar(y ~ x1 + x2, "deff", "n_i", weight_mat, data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 4: Auxiliary variable (X) contains NA values
  data_invalid_x <- databeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    betadeff_sar(y ~ x1 + x2, "deff", "n_i", weight_mat, data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 5: Iteration update is less than 3
  expect_error(
    betadeff_sar(y ~ x1 + x2, "deff", "n_i", weight_mat, databeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 6: Sample size (n_i) is less than or equal to deff
  data_invalid_n <- databeta
  data_invalid_n$n_i[1] <- 1.5
  data_invalid_n$deff[1] <- 2.0
  expect_error(
    betadeff_sar(y ~ x1 + x2, "deff", "n_i", weight_mat, data_invalid_n, plot = FALSE),
    "There is at least one sampled area where n_i <= deff. Effective sample size must be > 1"
  )

  # Case 7: Proximity matrix dimension mismatch
  wrong_W <- weight_mat[1:10, 1:10]
  expect_error(
    betadeff_sar(y ~ x1 + x2, "deff", "n_i", wrong_W, databeta, plot = FALSE),
    "Proximity matrix must be N x N"
  )

  # Case 8: Proximity matrix contains NA
  na_W <- weight_mat
  na_W[1, 2] <- NA
  expect_error(
    betadeff_sar(y ~ x1 + x2, "deff", "n_i", na_W, databeta, plot = FALSE),
    "Proximity matrix contains NA"
  )

  # Case 9: Formula without predictor
  expect_error(
    betadeff_sar(y ~ 1, "deff", "n_i", weight_mat, databeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
})
