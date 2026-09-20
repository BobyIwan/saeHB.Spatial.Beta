test_that("Unit Testing for beta_nonspatial: Execution and Output Format", {
  skip_on_cran()

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- beta_nonspatial(
      formula = y ~ x1 + x2,
      data = databeta,

      iter.mcmc = 100,
      burn.in = 50,
      chains = 2,
      quiet = TRUE
    )
    dev.off()
    expect_true(is.list(res_sampled))
    expect_equal(rownames(res_sampled$coefficient), c("beta[0]", "beta[1]", "beta[2]", "phi"))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- beta_nonspatial(
      formula = y ~ x1 + x2,
      data = databeta_na,

      iter.mcmc = 100,
      burn.in = 50,
      chains = 2,

      quiet = TRUE,
      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
    expect_equal(rownames(res_nonsampled$coefficient), c("beta[0]", "beta[1]", "beta[2]", "phi"))
  })
})

test_that("Unit Testing for beta_nonspatial: Error Handling", {
  # Case 3: Response variable (y) not between 0 and 1
  data_invalid_y <- databeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    beta_nonspatial(y ~ x1 + x2, data = data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 4: Auxiliary variable (X) contains NA values
  data_invalid_x <- databeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    beta_nonspatial(y ~ x1 + x2, data = data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 5: Formula without predictor or intercept
  expect_error(
    beta_nonspatial(y ~ 1, data = databeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
  expect_error(
    beta_nonspatial(y ~ x1 - 1, data = databeta, plot = FALSE),
    "Model must include an intercept"
  )

  # Case 6: MCMC and Prior hyperparameters validation errors
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, iter.update = 2, plot = FALSE), "The number of iteration updates must be at least 3")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, iter.mcmc = 50, burn.in = 100, plot = FALSE), "iter.mcmc must exceed burn.in")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, thin = 0, plot = FALSE), "thin must be >= 1")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, chains = 0, plot = FALSE), "chains must be >= 1")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, tau.v = -1, plot = FALSE), "tau.v must be positive")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, seed = 0, plot = FALSE), "seed must be positive")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, coef = c(1, 2), plot = FALSE), "coef must have length equal to")
  expect_error(beta_nonspatial(y ~ x1 + x2, data = databeta, var.coef = c(1, 1, -1), plot = FALSE), "All values in var.coef must be positive")

  # Case 7: n.sims validation errors
  expect_error(
    beta_nonspatial(y ~ x1 + x2, data = databeta, n.sims = 0, plot = FALSE),
    "n.sims must be >= 1"
  )
  expect_error(
    beta_nonspatial(y ~ x1 + x2, data = databeta, chains = 2, n.sims = 3, plot = FALSE),
    "n.sims cannot exceed the number of chains"
  )

  # Case 8: Effective sample size validation error for runjags
  expect_error(
    beta_nonspatial(y ~ x1 + x2, data = databeta, iter.mcmc = 100, burn.in = 99, thin = 2, plot = FALSE),
    "The effective number of samples"
  )
})
