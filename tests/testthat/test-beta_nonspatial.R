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
      quiet = TRUE
    )
    dev.off()
    expect_true(is.list(res_sampled))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- beta_nonspatial(
      formula = y ~ x1 + x2,
      data = databeta_na,

      iter.mcmc = 100,
      burn.in = 50,

      quiet = TRUE,
      plot = FALSE
    )
    expect_true(is.list(res_nonsampled))
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

  # Case 5: Iteration update is less than 3
  expect_error(
    beta_nonspatial(y ~ x1 + x2, data = databeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 6: Formula without predictor
  expect_error(
    beta_nonspatial(y ~ 1, data = databeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
})
