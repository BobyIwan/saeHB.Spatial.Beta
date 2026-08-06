test_that("Unit Testing for betadeff_lerouxcar: Execution and Output Format", {
  skip_on_cran()

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- betadeff_lerouxcar(
      formula = y ~ x1 + x2,
      deff = "deff",
      n_i = "n_i",
      proxmat = adjacency_mat,
      data = databeta,

      iter.mcmc = 100,
      burn.in = 50,
      quiet = TRUE,
      n.adapt = 2000
    )
    dev.off()
    expect_true(is.list(res_sampled))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- betadeff_lerouxcar(
      formula = y ~ x1 + x2,
      deff = "deff",
      n_i = "n_i",
      proxmat = adjacency_mat,
      data = databeta_na,

      iter.mcmc = 100,
      burn.in = 50,
      n.adapt = 2000,

      quiet = TRUE,
      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
  })
})

test_that("Unit Testing for betadeff_lerouxcar: Error Handling", {
  # Case 3: Response variable (y) not between 0 and 1
  data_invalid_y <- databeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 4: Auxiliary variable (X) contains NA values
  data_invalid_x <- databeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 5: Iteration update is less than 3
  expect_error(
    betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 6: Sample size (n_i) is less than or equal to deff
  data_invalid_n <- databeta
  data_invalid_n$n_i[1] <- 1.5
  data_invalid_n$deff[1] <- 2.0
  expect_error(
    betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_invalid_n, plot = FALSE),
    "There is at least one sampled area where n_i <= deff"
  )

  # Case 7: Proximity matrix dimension mismatch and NA
  wrong_W <- adjacency_mat[1:10, 1:10]
  expect_error(
    betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", wrong_W, databeta, plot = FALSE),
    "Proximity matrix must be N x N"
  )

  na_W <- adjacency_mat
  na_W[1, 2] <- NA
  expect_error(
    betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", na_W, databeta, plot = FALSE),
    "Proximity matrix contains NA"
  )

  # Case 8: Formula without predictor or intercept
  expect_error(
    betadeff_lerouxcar(y ~ 1, "deff", "n_i", adjacency_mat, databeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
  expect_error(
    betadeff_lerouxcar(y ~ x1 - 1, "deff", "n_i", adjacency_mat, databeta, plot = FALSE),
    "Model must include an intercept"
  )

  # Case 9: Matrix binary, symmetry, and diagonal validation errors
  non_binary_W <- adjacency_mat; non_binary_W[1, 2] <- 0.5; non_binary_W[2, 1] <- 0.5
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", non_binary_W, databeta, plot = FALSE), "proxmat must be a binary adjacency matrix containing only 0 and 1")

  asym_W <- adjacency_mat; asym_W[1, 2] <- 1; asym_W[2, 1] <- 0
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", asym_W, databeta, plot = FALSE), "proxmat must be symmetric")

  diag_W <- adjacency_mat; diag_W[1, 1] <- 1
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", diag_W, databeta, plot = FALSE), "Diagonal elements of proxmat must be zero")

  # Case 10: DEFF and n_i missing or negative value errors
  data_na_deff <- databeta; data_na_deff$deff[1] <- NA
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_na_deff, plot = FALSE), "Design effect contains NA values")

  data_neg_ni <- databeta; data_neg_ni$n_i[1] <- -1
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, data_neg_ni, plot = FALSE), "Sample sizes in sampled areas must be positive")

  # Case 11: MCMC and Prior hyperparameters validation errors
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, iter.mcmc = 50, burn.in = 100, plot = FALSE), "iter.mcmc must exceed burn.in")
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, thin = 0, plot = FALSE), "thin must be >= 1")
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, chains = 0, plot = FALSE), "chains must be >= 1")
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, tau.v = -1, plot = FALSE), "tau.v must be positive")
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, seed = 0, plot = FALSE), "seed must be positive")
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, coef = c(1, 2), plot = FALSE), "coef must have length equal to")
  expect_error(betadeff_lerouxcar(y ~ x1 + x2, "deff", "n_i", adjacency_mat, databeta, var.coef = c(1, 1, -1), plot = FALSE), "All values in var.coef must be positive")
})
