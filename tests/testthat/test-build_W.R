library(sf)

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

test_that("Unit Testing for build_W methods and formats", {
  suppressWarnings({
    # Case where methods execute successfully and output correct dimensions
    res_queen <- build_W(data = mock_poly, method = "contiguity", contiguity = "queen", output = "matrix")
    expect_true(is.matrix(res_queen))
    expect_equal(dim(res_queen), c(9, 9))

    res_knn <- build_W(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 2, output = "matrix")
    expect_true(is.matrix(res_knn))
    expect_equal(dim(res_knn), c(5, 5))

    res_kernel <- build_W(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = 50, output = "matrix")
    expect_true(is.matrix(res_kernel))
    expect_equal(dim(res_kernel), c(5, 5))

    # Case where format output variants defined
    expect_s3_class(build_W(data = mock_poly, method = "contiguity", output = "listw"), "listw")
    expect_type(build_W(data = mock_poly, method = "contiguity", output = "all"), "list")
    expect_true(inherits(build_W(data = mock_poly, method = "contiguity", output = "nb"), "nb"))
  })
})

test_that("Unit Testing for build_W variants and fallbacks", {
  suppressWarnings({
    # Case where contiguity is rook and bishop
    expect_equal(dim(build_W(data = mock_poly, method = "contiguity", contiguity = "rook", output = "matrix")), c(9, 9))
    expect_equal(dim(build_W(data = mock_poly, method = "contiguity", contiguity = "bishop", output = "matrix")), c(9, 9))

    # Case where fallback is knn and distance
    expect_true(sum(build_W(data = mock_iso, method = "contiguity", fallback = "knn", fallback_k = 1, output = "matrix")[3, ]) > 0)
    expect_true(sum(build_W(data = mock_iso, method = "contiguity", fallback = "distance", fallback_dmax = 2000, output = "matrix")[3, ]) > 0)

    # Case where distance is inverse and exponential
    expect_equal(dim(build_W(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", power = 2, output = "matrix")), c(5, 5))
    expect_equal(dim(build_W(data = mock_coords, coords = mock_coords, method = "distance", distance = "exponential", alpha = 0.5, output = "matrix")), c(5, 5))

    # Case where distance using dmax
    expect_equal(dim(build_W(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", dmax = 50, output = "matrix")), c(5, 5))
  })
})

test_that("Unit Testing for build_W kernel variants", {
  suppressWarnings({
    # Case where all kernel functions executed
    kernels <- c("uniform", "triangular", "epanechnikov", "quartic")
    for (k in kernels) {
      res_k <- build_W(data = mock_coords, coords = mock_coords, method = "kernel", kernel = k, bandwidth = 50, output = "matrix")
      expect_true(is.matrix(res_k))
      expect_equal(dim(res_k), c(5, 5))
    }
  })
})

test_that("Unit Testing for build_W warnings and errors", {
  # Case where coordinate mismatch warnings occur
  suppressWarnings(expect_warning(build_W(data = mock_coords, coords = mock_coords, method = "distance", lonlat = FALSE), "Coordinate data looks like"))
  suppressWarnings(expect_warning(build_W(data = planar_coords, coords = planar_coords, method = "distance", lonlat = TRUE), "Coordinate data does NOT look like"))

  # Case where zero policy is FALSE
  expect_error(suppressWarnings(build_W(data = mock_iso, method = "contiguity", fallback = "none", zero.policy = FALSE)))

  # Case where basic errors captured
  invalid_data <- data.frame(x = 1:5, y = 1:5)
  mock_na <- mock_coords; mock_na[1, 1] <- NA

  expect_error(build_W(data = invalid_data, method = "contiguity"))
  expect_error(build_W(data = invalid_data, method = "distance"))
  expect_error(build_W(data = mock_coords, coords = mock_coords, method = "kernel", bandwidth = NULL))
  expect_error(suppressWarnings(build_W(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 0)))
  expect_error(build_W(data = mock_na, coords = mock_na, method = "distance"))
})
