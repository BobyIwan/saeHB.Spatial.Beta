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
    # Case 1: Contiguity method (queen) executes successfully and outputs correct dimensions
    res_queen <- build_w(data = mock_poly, method = "contiguity", contiguity = "queen", output = "matrix")
    expect_true(is.matrix(res_queen))
    expect_equal(dim(res_queen), c(9, 9))

    # Case 2: Distance method (knn) executes successfully and outputs correct dimensions
    res_knn <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 2, output = "matrix")
    expect_true(is.matrix(res_knn))
    expect_equal(dim(res_knn), c(5, 5))

    # Case 3: Kernel method (gaussian) executes successfully and outputs correct dimensions
    res_kernel <- build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = 50, output = "matrix")
    expect_true(is.matrix(res_kernel))
    expect_equal(dim(res_kernel), c(5, 5))

    # Case 4: Format output variants (listw, all, nb) execute correctly
    expect_s3_class(build_w(data = mock_poly, method = "contiguity", output = "listw"), "listw")
    expect_type(build_w(data = mock_poly, method = "contiguity", output = "all"), "list")
    expect_true(inherits(build_w(data = mock_poly, method = "contiguity", output = "nb"), "nb"))
  })
})

test_that("Unit Testing for build_w variants and fallbacks", {
  suppressWarnings({
    # Case 5: Contiguity variants (rook and bishop) execute successfully
    expect_equal(dim(build_w(data = mock_poly, method = "contiguity", contiguity = "rook", output = "matrix")), c(9, 9))
    expect_equal(dim(build_w(data = mock_poly, method = "contiguity", contiguity = "bishop", output = "matrix")), c(9, 9))

    # Case 6: Fallback mechanisms (knn and distance) for isolated areas work correctly
    expect_true(sum(build_w(data = mock_iso, method = "contiguity", fallback = "knn", fallback_k = 1, output = "matrix")[3, ]) > 0)
    expect_true(sum(build_w(data = mock_iso, method = "contiguity", fallback = "distance", fallback_dmax = 2000, output = "matrix")[3, ]) > 0)

    # Case 7: Distance variants (inverse and exponential) with style = W and B
    expect_equal(dim(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", power = 2, output = "matrix", style = "W")), c(5, 5))
    res_exp_b <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "exponential", alpha = 0.5, output = "matrix", style = "B")
    expect_equal(dim(res_exp_b), c(5, 5))
    expect_true(all(res_exp_b %in% c(0, 1))) # Style B must be binary

    # Case 8: Distance method using dmax threshold executes successfully
    expect_equal(dim(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", dmax = 50, output = "matrix")), c(5, 5))
  })
})

test_that("Unit Testing for build_w: style='B' produces symmetric matrices", {
  suppressWarnings({
    # Case 9: Contiguity + fallback must produce symmetric W for style="B"
    W_contig_fallback_b <- build_w(data = mock_iso, method = "contiguity", fallback = "knn", fallback_k = 1, style = "B", output = "matrix")
    expect_true(isSymmetric(W_contig_fallback_b))

    # Case 10: distance + knn must produce symmetric W for style="B"
    W_knn_b <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 2, style = "B", output = "matrix")
    expect_true(isSymmetric(W_knn_b))

    # Case 11: inverse_distance without dmax must produce symmetric W for style="B"
    W_invdist_b <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", power = 2, style = "B", output = "matrix")
    expect_true(isSymmetric(W_invdist_b))

    # Case 12: exponential without dmax must produce symmetric W for style="B"
    W_exp_b <- build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "exponential", alpha = 0.5, style = "B", output = "matrix")
    expect_true(isSymmetric(W_exp_b))
  })
})

test_that("Unit Testing for build_w kernel variants", {
  suppressWarnings({
    # Case 13: All non-gaussian kernel functions execute successfully
    kernels <- c("uniform", "triangular", "epanechnikov", "quartic")
    for (k in kernels) {
      res_k <- build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = k, bandwidth = 50, output = "matrix")
      expect_true(is.matrix(res_k))
      expect_equal(dim(res_k), c(5, 5))
    }

    # Case 14: Kernel method with style = B produces binary matrix
    res_kernel_b <- build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = 50, output = "matrix", style = "B")
    expect_true(all(res_kernel_b %in% c(0, 1)))
  })
})

test_that("Unit Testing for build_w warnings and errors (Edge Cases)", {
  # Case 15: Warnings for coordinate projection mismatch
  suppressWarnings(expect_warning(build_w(data = mock_coords, coords = mock_coords, method = "distance", lonlat = FALSE), "Coordinate data looks like"))
  suppressWarnings(expect_warning(build_w(data = planar_coords, coords = planar_coords, method = "distance", lonlat = TRUE), "Coordinate data does NOT look like"))

  invalid_data <- data.frame(x = 1:5, y = 1:5)
  mock_na <- mock_coords; mock_na[1, 1] <- NA
  far_coords <- rbind(mock_coords, c(lon = 150, lat = 150)) # To create isolation

  # Case 16: Errors for invalid coordinates and dimensions
  expect_error(build_w(data = NULL, coords = c(1,2), method = "distance"), "coords must be numeric with 2 columns")
  expect_error(build_w(data = invalid_data, method = "contiguity"))
  expect_error(build_w(data = invalid_data, method = "distance"))
  expect_error(build_w(data = mock_na, coords = mock_na, method = "distance"), "Coordinates must be finite")

  # Case 17: Error when less than 2 areas are provided
  expect_error(build_w(data = mock_coords[1,,drop=FALSE], coords = mock_coords[1,,drop=FALSE], method = "distance"), "Need at least 2 areas")

  # Case 18: Error for missing fallback_dmax in contiguity
  expect_error(suppressWarnings(build_w(data = mock_iso, method = "contiguity", fallback = "distance", fallback_dmax = NULL)), "fallback_dmax must be provided")

  # Case 19: Errors for invalid distance parameters (k, power, alpha)
  expect_error(suppressWarnings(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "knn", k = 0)), "must be a single positive integer")
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "inverse_distance", power = -1), "must be a single positive numeric value")
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "distance", distance = "exponential", alpha = 0), "must be a single positive numeric value")

  # Case 20: Errors for invalid kernel parameters (bandwidth)
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "kernel", bandwidth = NULL), "must be a single positive numeric value")
  expect_error(build_w(data = mock_coords, coords = mock_coords, method = "kernel", kernel = "gaussian", bandwidth = -10), "must be a single positive numeric value")

  # Case 21: Errors when zero.policy = FALSE on isolated areas
  expect_error(suppressWarnings(build_w(data = mock_iso, method = "contiguity", fallback = "none", zero.policy = FALSE)))
  expect_error(suppressWarnings(build_w(data = far_coords, coords = far_coords, method = "distance", distance = "inverse_distance", dmax = 10, zero.policy = FALSE)), "Some areas have no neighbors with current dmax")
  expect_error(suppressWarnings(build_w(data = far_coords, coords = far_coords, method = "kernel", kernel = "gaussian", bandwidth = 10, zero.policy = FALSE)), "Some areas have no neighbors with current bandwidth")
})
