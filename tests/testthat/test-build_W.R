boundary <- sf::st_sfc(sf::st_point(c(0,0)), sf::st_point(c(3,3)))
grid_geom <- sf::st_make_grid(boundary, n = c(3, 3))
mock_poly <- sf::st_sf(id = 1:9, geometry = grid_geom)

mock_coords <- cbind(lon = c(106.8, 106.9, 107.0, 107.1, 107.2),
                     lat = c(-6.2, -6.1, -6.3, -6.2, -6.4))

planar_coords <- cbind(x = c(500000, 501000, 502000), y = c(9000000, 9001000, 9002000))

p1 <- sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
p2 <- sf::st_polygon(list(rbind(c(1,0), c(2,0), c(2,1), c(1,1), c(1,0))))
p3 <- sf::st_polygon(list(rbind(c(10,10), c(11,10), c(11,11), c(10,11), c(10,10))))
mock_iso <- sf::st_sf(id=1:3, geometry=sf::st_sfc(p1,p2,p3))

test_that("Unit Testing for build_w methods and formats", {
  suppressWarnings({
    # Case where methods execute successfully and output correct dimensions
    res_queen <- build_w(data = mock_poly, method = "contiguity", contiguity = "queen", output = "matrix")
    expect_true(is.matrix(res_queen))
    expect_equal(dim(res_queen), c(9, 9))

    res_knn <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 2, output = "matrix")
    expect_true(is.matrix(res_knn))
    expect_equal(dim(res_knn), c(5, 5))

    res_kernel <- build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = 50, output = "matrix")
    expect_true(is.matrix(res_kernel))
    expect_equal(dim(res_kernel), c(5, 5))

    # Case where format output variants defined
    expect_s3_class(build_w(data = mock_poly, method = "contiguity", output = "listw"), "listw")
    expect_type(build_w(data = mock_poly, method = "contiguity", output = "all"), "list")
    expect_true(inherits(build_w(data = mock_poly, method = "contiguity", output = "nb"), "nb"))
  })
})

test_that("Unit Testing for build_w variants and fallbacks", {
  suppressWarnings({
    # Case where contiguity is rook and bishop
    expect_equal(dim(build_w(data = mock_poly, method = "contiguity", contiguity = "rook", output = "matrix")), c(9, 9))
    expect_equal(dim(build_w(data = mock_poly, method = "contiguity", contiguity = "bishop", output = "matrix")), c(9, 9))

    # Case where fallback is knn and distance
    expect_true(sum(build_w(data = mock_iso, method = "contiguity", fallback = "knn", fallback_k = 1, output = "matrix")[3, ]) > 0)
    expect_true(sum(build_w(data = mock_iso, method = "contiguity", fallback = "distance", fallback_dmax = 2000, output = "matrix")[3, ]) > 0)

    # Case where distance is inverse and exponential with style = W and B
    expect_equal(dim(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", power = 2, output = "matrix", style = "W")), c(5, 5))
    res_exp_b <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "exponential", alpha = 0.5, output = "matrix", style = "B")
    expect_equal(dim(res_exp_b), c(5, 5))
    expect_true(all(res_exp_b %in% c(0, 1))) # Style B must be binary

    # Case where distance using dmax
    expect_equal(dim(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", dmax = 50, output = "matrix")), c(5, 5))
  })
})

test_that("Unit Testing for build_w kernel variants", {
  suppressWarnings({
    # Case where all kernel functions executed
    kernels <- c("uniform", "triangular", "epanechnikov", "quartic")
    for (k in kernels) {
      res_k <- build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = k, bandwidth = 50, output = "matrix")
      expect_true(is.matrix(res_k))
      expect_equal(dim(res_k), c(5, 5))
    }

    # Test kernel style = B
    res_kernel_b <- build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = 50, output = "matrix", style = "B")
    expect_true(all(res_kernel_b %in% c(0, 1))) # Style B must be binary
  })
})

test_that("Unit Testing for build_w warnings and errors (Edge Cases)", {
  # Warning cases
  suppressWarnings(expect_warning(build_w(data = mock_coords, coords = mock_coords, method = "distance", lonlat = FALSE), "Coordinate data looks like"))
  suppressWarnings(expect_warning(build_w(data = planar_coords, coords = planar_coords, method = "distance", lonlat = TRUE), "Coordinate data does NOT look like"))

  invalid_data <- data.frame(x = 1:5, y = 1:5)
  mock_na <- mock_coords; mock_na[1, 1] <- NA
  far_coords <- rbind(mock_coords, c(lon = 150, lat = 150)) # To create isolation

  # Error: Coords validation
  expect_error(build_w(data = NULL, coords = c(1,2), method = "distance"), "coords must be numeric with 2 columns")
  expect_error(build_w(data = invalid_data, method = "contiguity"))
  expect_error(build_w(data = invalid_data, method = "distance"))
  expect_error(build_w(data = mock_na, coords = mock_na, method = "distance"), "Some areas have empty geometries")

  # Error: Less than 2 areas
  expect_error(build_w(data = mock_coords[1,,drop=FALSE], coords = mock_coords[1,,drop=FALSE], method = "distance"), "Need at least 2 areas")

  # Error: Contiguity fallback without dmax (DIBUNGKAM DI SINI)
  expect_error(suppressWarnings(build_w(data = mock_iso, method = "contiguity", fallback = "distance", fallback_dmax = NULL)), "fallback_dmax must be provided")

  # Error: Distance invalid parameters
  expect_error(suppressWarnings(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 0)), "k must be >= 1")
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", power = -1), "power must be > 0")
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "exponential", alpha = 0), "alpha must be > 0")

  # Error: Kernel invalid parameters
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "kernel", bandwidth = NULL), "bandwidth must be provided")
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = -10), "bandwidth must be provided and > 0")

  # Error: Zero policy = FALSE on isolated areas
  expect_error(suppressWarnings(build_w(data = mock_iso, method = "contiguity", fallback = "none", zero.policy = FALSE)))
  expect_error(suppressWarnings(build_w(data = far_coords, coords = far_coords, method = "distance", distance = "inverse_distance", dmax = 10, zero.policy = FALSE)), "Some areas have no neighbors with current dmax")
  expect_error(suppressWarnings(build_w(data = far_coords, coords = far_coords, method = "kernel", kernel = "gaussian", bandwidth = 10, zero.policy = FALSE)), "Some areas have no neighbors with current bandwidth")
})
