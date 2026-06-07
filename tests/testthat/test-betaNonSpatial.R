test_that("Unit Testing for betaNonSpatial: Execution and Output Format", {

  suppressWarnings({
    # Case 1: Fully sampled data
    pdf(file = NULL)
    res_sampled <- betaNonSpatial(
      formula = y ~ x1 + x2,
      data = dataBeta,
      plot = TRUE,
      keep.fit = TRUE
    )
    dev.off()

    expect_true(is.list(res_sampled))
    expect_true(all(c("Est", "refVar", "randeff", "coefficient", "fit") %in% names(res_sampled)))
    expect_equal(nrow(res_sampled$Est), nrow(dataBeta))

    # Case 2: Data with non-sampled areas (NA) executes successfully
    res_nonsampled <- betaNonSpatial(
      formula = y ~ x1 + x2,
      data = dataBeta_NA,
      plot = FALSE
    )

    expect_true(is.list(res_nonsampled))
    expect_equal(nrow(res_nonsampled$Est), nrow(dataBeta_NA))

    # Case 3: Execution works with explicitly defined coef and var.coef
    res_coef <- betaNonSpatial(
      formula = y ~ x1 + x2,
      data = dataBeta,
      coef = c(0, 0, 0),
      var.coef = c(1, 1, 1),
      plot = FALSE
    )
    expect_true(is.list(res_coef))
  })
})

test_that("Unit Testing for betaNonSpatial: Error Handling", {
  # Case 4: Response variable (y) not between 0 and 1
  data_invalid_y <- dataBeta
  data_invalid_y$y[5] <- 1.5
  expect_error(
    betaNonSpatial(y ~ x1 + x2, data = data_invalid_y, plot = FALSE),
    "Response variable must satisfy 0 < y < 1"
  )

  # Case 5: Auxiliary variable (X) contains NA values
  data_invalid_x <- dataBeta
  data_invalid_x$x1[10] <- NA
  expect_error(
    betaNonSpatial(y ~ x1 + x2, data = data_invalid_x, plot = FALSE),
    "Auxiliary variables contain NA values"
  )

  # Case 6: Iteration update is less than 3
  expect_error(
    betaNonSpatial(y ~ x1 + x2, data = dataBeta, iter.update = 2, plot = FALSE),
    "The number of iteration updates must be at least 3"
  )

  # Case 7: Formula without predictor
  expect_error(
    betaNonSpatial(y ~ 1, data = dataBeta, plot = FALSE),
    "Formula must include response and at least 1 predictor"
  )
})
