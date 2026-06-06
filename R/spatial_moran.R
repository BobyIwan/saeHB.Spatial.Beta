#' @title Moran's I Test for Spatial Autocorrelation
#'
#' @description
#' \itemize{
#'   \item {This function provides a convenient wrapper to perform Moran's I test for spatial autocorrelation on a numeric vector.}
#'   \item {It seamlessly handles missing values (NA) by subsetting both the numeric vector and the spatial weights list simultaneously.}
#'   \item {It supports both analytical (randomization) and Monte Carlo permutation approaches.}
#' }
#'
#' @param x A numeric vector of the variable of interest (e.g., residuals, random effects, or raw data).
#' @param listw A \code{listw} object containing spatial weights created by \code{build_W} or \code{spdep}.
#' @param alternative A character string specifying the alternative hypothesis. Must be one of \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.
#' @param mc Logical; if \code{TRUE}, performs Moran's I test using Monte Carlo permutations. Default is \code{FALSE} (analytical approach).
#' @param nsim An integer specifying the number of permutations if \code{mc = TRUE}. Default is \code{999}.
#' @param zero.policy Logical; if \code{TRUE}, allows areas with no neighbors (isolates) to be included in the calculation. Default is \code{TRUE}.
#' @param na.rm Logical; if \code{TRUE}, missing values in \code{x} are removed, and the corresponding rows/columns in the spatial weights are automatically subsetted. Default is \code{TRUE}.
#'
#' @return A list with class \code{htest} containing the following components:
#' \itemize{
#'   \item \code{statistic}: The value of the standard deviate of Moran's I.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{estimate}: The value of the observed Moran's I, its expectation, and variance.
#'   \item \code{method}: A character string indicating the type of test performed.
#'   \item \code{data.name}: A character string giving the name(s) of the data.
#' }
#'
#' @import spdep
#'
#' @export spatial_moran
spatial_moran <- function(x,
                          listw,
                          alternative = c("greater", "less", "two.sided"),
                          mc = FALSE,
                          nsim = 999,
                          zero.policy = TRUE,
                          na.rm = TRUE) {

  alternative <- match.arg(alternative)

  var_name <- deparse(substitute(x))

  if (!is.numeric(x)) stop("Argument 'x' must be a numeric vector.")
  if (!inherits(listw, "listw")) {
    stop("Argument 'listw' must be an object of class 'listw'.")
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
    res$method <- "Moran's I test under randomization"
  }

  res$data.name <- var_name

  return(res)
}
