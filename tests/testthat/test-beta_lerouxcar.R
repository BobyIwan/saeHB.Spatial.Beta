test_that("Unit Testing for beta_lerouxcar: Execution and Output Format", {
  skip_on_cran()

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- beta_lerouxcar(
      formula = y ~ x1 + x2,
      proxmat = adjacency_mat,
      data = databeta,

      iter.mcmc = 100,
      burn.in = 50,
      quiet = TRUE,
      n.adapt = 2500
    )
    dev.off()
    expect_true(is.list(res_sampled))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- beta_lerouxcar(
      formula = y ~ x1 + x2,
      proxmat = adjacency_mat,
      data = databeta_na,

      iter.mcmc = 100,
      burn.in = 50,
      quiet = TRUE,
      n.adapt = 2500,

      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
  })
})

test_that("Unit Testing for beta_lerouxcar: Error Handling", {
  # Case 3: Response variable (y) not between 0 and 1
  data_invalid_y <- databeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    beta_lerouxcar(y ~ x1 + x2, proxmat = adjacency_mat, data = data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 4: Auxiliary variable (X) contains NA values
  data_invalid_x <- databeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    beta_lerouxcar(y ~ x1 + x2, proxmat = adjacency_mat, data = data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 5: Iteration update is less than 3
  expect_error(
    beta_lerouxcar(y ~ x1 + x2, proxmat = adjacency_mat, data = databeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 6: Proximity matrix dimension mismatch
  wrong_W <- adjacency_mat[1:10, 1:10]
  expect_error(
    beta_lerouxcar(y ~ x1 + x2, proxmat = wrong_W, data = databeta, plot = FALSE),
    "Proximity matrix must be N x N"
  )

  # Case 7: Proximity matrix contains NA
  na_W <- adjacency_mat
  na_W[1, 2] <- NA
  expect_error(
    beta_lerouxcar(y ~ x1 + x2, proxmat = na_W, data = databeta, plot = FALSE),
    "Proximity matrix contains NA"
  )

  # Case 8: Formula without predictor
  expect_error(
    beta_lerouxcar(y ~ 1, proxmat = adjacency_mat, data = databeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
})
