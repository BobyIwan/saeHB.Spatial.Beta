#' @title Moran's I Test for Spatial Autocorrelation
#'
#' @description This function performs Moran's I test for detecting global spatial autocorrelation. It provides a convenient wrapper around \code{spdep::moran.test()} and \code{spdep::moran.mc()}, with automatic handling of missing values (NA) by seamlessly subsetting both the response vector and the spatial weights object simultaneously.
#'
#' @details
#' This function supports two approaches to testing the significance of Moran's I:
#'
#' \strong{1. Analytical Approach (Randomization - Default)}
#' \cr
#' When \code{mc = FALSE}, the function uses the analytical approach (specifically, the assumption of randomization). It computes the theoretical expectation and variance of Moran's I under the null hypothesis of no spatial autocorrelation. This approach relies on an asymptotic approximation of the sampling distribution of Moran's I.
#' \cr
#' \emph{When to use:} Use this approach when your dataset is relatively large and follows standard statistical assumptions. It is computationally fast and provides reliable asymptotic p-values for large \eqn{N}.
#'
#' \strong{2. Monte Carlo Permutation Approach (\code{mc = TRUE})}
#' \cr
#' When \code{mc = TRUE}, the function calculates the p-value empirically. It randomly permutes (shuffles) the observed values \code{x} across the spatial units \code{nsim} times. Because it computes the p-value empirically without relying on asymptotic theory, the Monte Carlo permutation approach is particularly useful for small datasets or when the assumptions of the analytical test may not hold.
#'
#' @param x A numeric vector of the variable of interest (e.g., residuals, random effects, or raw data).
#' @param listw A \code{listw} object containing spatial weights, typically created by \code{build_w()}.
#' @param alternative A character string specifying the alternative hypothesis. Must be one of \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.
#' @param mc Logical; if \code{TRUE}, performs Moran's I test using Monte Carlo permutations. Default is \code{FALSE} (analytical approach).
#' @param nsim An integer specifying the number of permutations if \code{mc = TRUE}. Default is \code{999}.
#' @param zero.policy Logical; if \code{TRUE}, allows areas with no neighbors (isolates) to be included in the calculation. Default is \code{TRUE}.
#' @param na.rm Logical; if \code{TRUE}, missing values in \code{x} are removed, and the corresponding rows/columns in the spatial weights are automatically subsetted. Default is \code{TRUE}.
#'
#' @return An object of class \code{htest} (if \code{mc = FALSE}) or \code{mc.sim} (if \code{mc = TRUE}). The returned components depend on the selected test, following the corresponding \code{spdep} implementation. Common components include:
#' \itemize{
#'   \item \code{statistic}: The value of the standard deviate of Moran's I.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{method}: A character string indicating the type of test performed.
#'   \item \code{data.name}: A character string giving the name(s) of the data.
#' }
#'
#' @examples
#' library(sf)
#'
#' # 1. Prepare dummy data
#' bbox <- st_bbox(c(xmin = 0, ymin = 0, xmax = 3, ymax = 3))
#' grid <- st_make_grid(bbox, n = c(3, 3))
#' grid_sf <- st_sf(id = 1:9, geometry = grid)
#' set.seed(123)
#' grid_sf$y <- rnorm(9)
#'
#' # 2. Build spatial weights using the package's native function
#' W_obj <- build_w(
#'   data = grid_sf,
#'   method = "contiguity",
#'   contiguity = "queen",
#'   output = "all"
#' )
#'
#' # 3. Perform Moran's I test (Analytical approach)
#' moran_test(x = grid_sf$y, listw = W_obj$listw)
#'
#' # 4. Perform Moran's I test (Monte Carlo permutation approach)
#' moran_test(x = grid_sf$y, listw = W_obj$listw, mc = TRUE, nsim = 99)
#'
#' # 5. Handling Missing Values automatically
#' y_with_na <- grid_sf$y
#' y_with_na[c(2, 5)] <- NA
#' moran_test(x = y_with_na, listw = W_obj$listw, na.rm = TRUE)
#'
#' @import spdep
#'
#' @export moran_test
moran_test <- function(x,
                       listw,
                       alternative = c("greater", "less", "two.sided"),
                       mc = FALSE,
                       nsim = 999,
                       zero.policy = TRUE,
                       na.rm = TRUE) {

  alternative <- match.arg(alternative)
  var_name <- deparse(substitute(x))

  if (!is.numeric(x)) stop("Argument 'x' must be a numeric vector.")
  if (!inherits(listw, "listw")) stop("Argument 'listw' must be an object of class 'listw'.")

  # Validation: Check if length of x matches the spatial weights dimensions
  if (length(x) != length(listw$neighbours)) {
    stop(
      sprintf(
        "Length of x (%d) does not match the number of spatial units in listw (%d).",
        length(x),
        length(listw$neighbours)
      )
    )
  }

  if (any(is.na(x))) {
    if (na.rm) {
      valid_idx <- !is.na(x)
      n_valid <- sum(valid_idx)

      if (n_valid < 3) {
        stop("Not enough valid data points (non-NA) to compute Moran's I. Minimum is 3.")
      }

      n_missing <- length(x) - n_valid
      x_clean <- x[valid_idx]

      listw <- tryCatch({
        spdep::subset.listw(listw, subset = valid_idx, zero.policy = zero.policy)
      }, error = function(e) {
        W_mat <- spdep::listw2mat(listw)
        W_mat_sub <- W_mat[valid_idx, valid_idx, drop = FALSE]
        style <- if (!is.null(listw$style)) listw$style else "W"
        spdep::mat2listw(W_mat_sub, style = style, zero.policy = zero.policy)
      })

      message(sprintf("Info: %d missing values detected and removed. Spatial weights subsetted accordingly.", n_missing))
      var_name <- paste0(var_name, " (", n_missing, " NA removed)")
    } else {
      stop("Missing values (NA) detected in 'x'. Set na.rm = TRUE to handle them automatically.")
    }
  } else {
    x_clean <- x
  }

  if (mc) {
    res <- spdep::moran.mc(x = x_clean,
                           listw = listw,
                           nsim = nsim,
                           alternative = alternative,
                           zero.policy = zero.policy)
    res$method <- "Moran's I test under Monte Carlo permutation"
  } else {
    res <- spdep::moran.test(x = x_clean,
                             listw = listw,
                             alternative = alternative,
                             zero.policy = zero.policy)
    res$method <- "Moran's I test under randomization (analytical)"
  }

  res$data.name <- var_name

  return(res)
}
