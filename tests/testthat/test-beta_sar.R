test_that("Unit Testing for beta_sar: Execution and Output Format", {
  skip_on_cran()

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- beta_sar(
      formula = y ~ x1 + x2,
      proxmat = weight_mat,
      data = databeta,

      iter.mcmc = 100,
      burn.in = 50,
      quiet = TRUE
    )
    dev.off()
    expect_true(is.list(res_sampled))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- beta_sar(
      formula = y ~ x1 + x2,
      proxmat = weight_mat,
      data = databeta_na,

      iter.mcmc = 100,
      burn.in = 50,

      quiet = TRUE,
      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
  })
})

test_that("Unit Testing for beta_sar: Error Handling", {
  # Case 3: Response variable (y) not between 0 and 1
  data_invalid_y <- databeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 4: Auxiliary variable (X) contains NA values
  data_invalid_x <- databeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 5: Formula without predictor or intercept
  expect_error(
    beta_sar(y ~ 1, proxmat = weight_mat, data = databeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
  expect_error(
    beta_sar(y ~ x1 - 1, proxmat = weight_mat, data = databeta, plot = FALSE),
    "Model must include an intercept"
  )

  # Case 6: Proximity matrix dimension mismatch and NA
  wrong_W <- weight_mat[1:10, 1:10]
  expect_error(
    beta_sar(y ~ x1 + x2, proxmat = wrong_W, data = databeta, plot = FALSE),
    "Proximity matrix must be N x N"
  )

  na_W <- weight_mat
  na_W[1, 2] <- NA
  expect_error(
    beta_sar(y ~ x1 + x2, proxmat = na_W, data = databeta, plot = FALSE),
    "Proximity matrix contains NA"
  )

  # Case 7: Matrix row-standardization and diagonal validation errors
  unstd_W <- weight_mat; unstd_W[1, 2] <- 0.99
  expect_error(
    beta_sar(y ~ x1 + x2, proxmat = unstd_W, data = databeta, plot = FALSE),
    "Each row of proxmat must sum to 1"
  )

  diag_W <- weight_mat; diag_W[1, 1] <- 0.5; diag_W[1, 2] <- diag_W[1, 2] - 0.5
  expect_error(
    beta_sar(y ~ x1 + x2, proxmat = diag_W, data = databeta, plot = FALSE),
    "Diagonal elements of proxmat must be zero"
  )

  # Case 8: MCMC and Prior hyperparameters validation errors
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, iter.update = 2, plot = FALSE), "The number of iteration updates must be at least 3")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, iter.mcmc = 50, burn.in = 100, plot = FALSE), "iter.mcmc must exceed burn.in")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, thin = 0, plot = FALSE), "thin must be >= 1")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, chains = 0, plot = FALSE), "chains must be >= 1")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, tau.u = -1, plot = FALSE), "tau.u must be positive")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, seed = 0, plot = FALSE), "seed must be positive")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, coef = c(1, 2), plot = FALSE), "coef must have length equal to")
  expect_error(beta_sar(y ~ x1 + x2, proxmat = weight_mat, data = databeta, var.coef = c(1, 1, -1), plot = FALSE), "All values in var.coef must be positive")
})
