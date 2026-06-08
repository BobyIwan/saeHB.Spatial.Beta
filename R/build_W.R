#' @title Build Spatial Weights Matrix
#'
#' @description
#' \itemize{
#'   \item {This function constructs spatial weights matrices (\eqn{W}) for spatial modeling.}
#'   \item {It supports various methods including Contiguity (Queen, Rook, Bishop), Distance-based (K-Nearest Neighbors, Inverse Distance, Exponential), and Kernel-based weights.}
#'   \item {It also provides a robust fallback mechanism to automatically connect isolated areas (islands) when using contiguity methods.}
#' }
#'
#' @param data An \code{sf} object (polygons/points) or a standard data frame. If a standard data frame is provided, \code{coords} must be specified.
#' @param coords An \eqn{N \times 2} matrix of coordinates. Required if \code{data} is not an \code{sf} object.
#' @param method A string indicating the spatial weight construction method. Options are \code{"contiguity"}, \code{"distance"}, or \code{"kernel"}.
#' @param contiguity A string indicating the contiguity type. Options are \code{"queen"}, \code{"rook"}, or \code{"bishop"}.
#' @param distance A string indicating the distance-based type. Options are \code{"knn"}, \code{"inverse_distance"}, or \code{"exponential"}.
#' @param k An integer specifying the number of nearest neighbors for KNN methods. Default is \code{2}.
#' @param dmax A numeric specifying the maximum distance threshold for distance-based neighbors.
#' @param power A numeric specifying the decay power for inverse distance weights. Default is \code{1}.
#' @param alpha A numeric specifying the decay parameter for exponential distance weights. Default is \code{1}.
#' @param epsilon A small numeric value to prevent division by zero in inverse distance calculation. Default is \code{1e-12}.
#' @param kernel A string indicating the type of spatial kernel. Options are \code{"uniform"}, \code{"gaussian"}, \code{"triangular"}, \code{"epanechnikov"}, or \code{"quartic"}.
#' @param bandwidth A numeric specifying the bandwidth (\eqn{h}) for kernel weights. Required if \code{method = "kernel"}.
#' @param lonlat Logical; if \code{TRUE}, coordinates are treated as Longitude/Latitude and great-circle spherical distances are calculated. Default is \code{TRUE}.
#' @param style A character string specifying the spatial weights coding scheme (\code{"W"} for row-standardized or \code{"B"} for binary). Default is \code{"W"}.
#' @param zero.policy Logical; if \code{TRUE}, areas with no neighbors are allowed to have zero-weight rows. Default is \code{TRUE}.
#' @param fallback A string indicating the fallback method for isolated areas (without neighbors) when using contiguity. Options are \code{"knn"}, \code{"distance"}, or \code{"none"}. Default is \code{"knn"}.
#' @param fallback_k An integer specifying the number of neighbors for the fallback method. Default is \code{2}.
#' @param fallback_dmax A numeric specifying the maximum distance for the fallback method.
#' @param output A string specifying the format of the output. Options are \code{"all"} (returns a comprehensive list), \code{"matrix"}, \code{"listw"}, or \code{"nb"}. Default is \code{"all"}.
#'
#' @return Depending on the \code{output} argument, this function returns:
#' \itemize{
#'   \item \code{"matrix"}: An \eqn{N \times N} spatial weights matrix.
#'   \item \code{"listw"}: A \code{listw} object compatible with \code{spdep} functions.
#'   \item \code{"nb"}: An \code{nb} (neighborhood) object.
#'   \item \code{"all"}: A list containing \code{W} (matrix), \code{listw}, \code{nb}, \code{info} (method details), and \code{diag} (diagnostic metrics for isolates and fallback).
#' }
#'
#' @examples
#' # Generate random Longitude and Latitude coordinates for 10 areas
#' # (e.g., somewhere roughly in Indonesia)
#' set.seed(123)
#' lon <- runif(10, min = 100, max = 140) # Longitude
#' lat <- runif(10, min = -10, max = 10)  # Latitude
#' coords <- cbind(lon, lat)
#'
#' # Build KNN distance-based weights (k = 2) using spherical distance
#' W_knn <- build_W(
#'   data = NULL,
#'   coords = coords,
#'   method = "distance",
#'   distance = "knn",
#'   k = 2,
#'   lonlat = TRUE,
#'   output = "matrix"
#' )
#'
#' # Build Gaussian Kernel weights using 500 km bandwidth
#' W_kernel <- build_W(
#'   data = NULL,
#'   coords = coords,
#'   method = "kernel",
#'   kernel = "gaussian",
#'   bandwidth = 500,
#'   lonlat = TRUE,
#'   output = "matrix"
#' )
#'
#' @import sf
#' @import spdep
#' @importFrom stats na.omit
#'
#' @export build_W
build_W <- function(
    data,
    coords = NULL,
    method = c("contiguity", "distance", "kernel"),
    contiguity = c("queen", "rook", "bishop"),
    distance = c("knn", "inverse_distance", "exponential"),
    k = 2,
    dmax = NULL,
    power = 1,
    alpha = 1,
    epsilon = 1e-12,
    kernel = c("uniform", "gaussian", "triangular", "epanechnikov", "quartic"),
    bandwidth = NULL,
    lonlat = TRUE,
    style = "W",
    zero.policy = TRUE,
    fallback = c("knn", "distance", "none"),
    fallback_k = 2,
    fallback_dmax = NULL,
    output = c("all", "matrix", "listw", "nb")
) {

  method     <- match.arg(method)
  contiguity <- match.arg(contiguity)
  distance   <- match.arg(distance)
  kernel     <- match.arg(kernel)
  fallback   <- match.arg(fallback)
  output     <- match.arg(output)

  get_coords <- function(data, coords) {
    if (!is.null(coords)) {
      cm <- as.matrix(coords)
      if (!is.numeric(cm) || ncol(cm) != 2) stop("coords must be numeric with 2 columns.")
      return(cm)
    }
    if (inherits(data, "sf")) {
      suppressWarnings(cc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(data))))
      if (ncol(cc) < 2) stop("Failed to extract centroid coordinates from sf.")
      return(as.matrix(cc[, 1:2, drop = FALSE]))
    }
    stop("Provide `coords` when `data` is not an sf object.")
  }

  is_sf <- inherits(data, "sf")
  if (method == "contiguity" && !is_sf) {
    stop("method='contiguity' requires `data` as sf polygons.")
  }

  coords_mat <- get_coords(data, coords)

  if (any(is.na(coords_mat))) {
    stop("Some areas have empty geometries or NA coordinates (Ghost Regions). Please clean your spatial data first.")
  }

  n <- nrow(coords_mat)
  if (n < 2) stop("Need at least 2 areas to build W.")

  k <- min(k, n - 1)
  fallback_k <- min(fallback_k, n - 1)

  is_looks_like_degree <- all(coords_mat[,1] >= -180 & coords_mat[,1] <= 180) &&
    all(coords_mat[,2] >= -90  & coords_mat[,2] <= 90)
  if (is_looks_like_degree && !lonlat) {
    warning("Coordinate data looks like Latitude/Longitude (degrees), but `lonlat = FALSE` is set. Calculating Euclidean distance on degree data can produce flawed spatial weights. Consider setting `lonlat = TRUE`.")
  }
  if (!is_looks_like_degree && lonlat) {
    warning("Coordinate data does NOT look like Latitude/Longitude (values exceed -180/180 or -90/90), but `lonlat = TRUE` is set. Calculating spherical distance on planar/projected data can produce flawed spatial weights. Consider setting `lonlat = FALSE`.")
  }

  nb <- NULL
  lw <- NULL
  W  <- NULL
  diag <- list()

  if (method == "contiguity") {

    if (contiguity == "queen") {
      nb_main <- spdep::poly2nb(data, queen = TRUE)
    } else if (contiguity == "rook") {
      nb_main <- spdep::poly2nb(data, queen = FALSE)
    } else if (contiguity == "bishop") {
      nb_q <- spdep::poly2nb(data, queen = TRUE)
      nb_r <- spdep::poly2nb(data, queen = FALSE)
      nb_main <- spdep::diffnb(nb_q, nb_r)
    }

    isolates <- which(spdep::card(nb_main) == 0)
    diag$isolates_before <- isolates
    diag$n_isolates_before <- length(isolates)

    if (length(isolates) > 0 && fallback != "none") {
      if (fallback == "knn") {
        knn <- spdep::knearneigh(coords_mat, k = fallback_k, longlat = lonlat)
        nb_fb <- spdep::knn2nb(knn)
      } else {
        if (is.null(fallback_dmax)) stop("fallback_dmax must be provided when fallback='distance'.")
        nb_fb <- spdep::dnearneigh(coords_mat, 0, fallback_dmax, longlat = lonlat)
      }

      for (i in isolates) {
        nb_main[[i]] <- nb_fb[[i]]
      }
    }

    nb <- nb_main
    diag$isolates_after <- which(spdep::card(nb) == 0)
    diag$n_isolates_after <- length(diag$isolates_after)

    lw <- spdep::nb2listw(nb, style = style, zero.policy = zero.policy)
    W  <- spdep::listw2mat(lw)
  }

  if (method == "distance") {

    if (distance == "knn") {
      if (!is.numeric(k) || k < 1) stop("k must be >= 1 for knn.")
      knn <- spdep::knearneigh(coords_mat, k = k, longlat = lonlat)
      nb  <- spdep::knn2nb(knn)

      diag$isolates_after <- which(spdep::card(nb) == 0)
      diag$n_isolates_after <- length(diag$isolates_after)

      lw <- spdep::nb2listw(nb, style = style, zero.policy = zero.policy)
      W  <- spdep::listw2mat(lw)
    }

    if (distance %in% c("inverse_distance", "exponential")) {
      if (distance == "inverse_distance" && (!is.numeric(power) || power <= 0)) stop("power must be > 0.")
      if (distance == "exponential" && (!is.numeric(alpha) || alpha <= 0)) stop("alpha must be > 0.")

      if (!is.null(dmax)) {
        nb <- spdep::dnearneigh(coords_mat, 0, dmax, longlat = lonlat)
        if (any(spdep::card(nb) == 0) && !zero.policy) {
          stop("Some areas have no neighbors with current dmax. Increase dmax or use k-based distance.")
        }
      } else {
        if (!is.numeric(k) || k < 1) stop("k must be >= 1 when dmax is NULL.")
        knn <- spdep::knearneigh(coords_mat, k = k, longlat = lonlat)
        nb  <- spdep::knn2nb(knn)
      }

      dist_list <- spdep::nbdists(nb, coords_mat, longlat = lonlat)

      if (distance == "inverse_distance") {
        glist <- lapply(dist_list, function(d) 1 / (pmax(d, epsilon)^power))
        diag$invdist_power <- power
      } else if (distance == "exponential") {
        glist <- lapply(dist_list, function(d) exp(-alpha * d))
        diag$exp_alpha <- alpha
      }

      lw <- spdep::nb2listw(nb, glist = glist, style = style, zero.policy = zero.policy)
      W  <- spdep::listw2mat(lw)

      diag$used_dmax <- !is.null(dmax)
      diag$dmax <- dmax
      diag$k <- k
    }
  }

  if (method == "kernel") {
    if (is.null(bandwidth) || !is.numeric(bandwidth) || bandwidth <= 0) {
      stop("bandwidth must be provided and > 0 for kernel weights.")
    }

    pts <- sf::st_as_sf(as.data.frame(coords_mat), coords = 1:2,
                        crs = ifelse(lonlat, 4326, NA))
    dist_mat <- as.numeric(sf::st_distance(pts))
    dist_mat <- matrix(dist_mat, n, n)
    if (lonlat == TRUE) {
      dist_mat <- dist_mat / 1000
    }

    Z <- dist_mat / bandwidth
    W <- matrix(0, n, n)

    if (kernel == "uniform") {
      W[dist_mat > 0 & dist_mat < bandwidth] <- 0.5
    }
    else if (kernel == "triangular") {
      idx <- which(dist_mat > 0 & dist_mat < bandwidth)
      W[idx] <- 1 - abs(Z[idx])
    }
    else if (kernel == "epanechnikov") {
      idx <- which(dist_mat > 0 & dist_mat < bandwidth)
      W[idx] <- 0.75 * (1 - Z[idx]^2)
    }
    else if (kernel == "quartic") {
      idx <- which(dist_mat > 0 & dist_mat < bandwidth)
      W[idx] <- (15/16) * (1 - Z[idx]^2)^2
    }
    else if (kernel == "gaussian") {
      W <- (1 / sqrt(2 * pi)) * exp(-(Z^2) / 2)
      diag(W) <- 0 # Tidak ada self-neighbor
    }

    if (any(rowSums(W) == 0) && !zero.policy) {
      stop("Some areas have no neighbors with current bandwidth. Increase bandwidth (h) or set zero.policy = TRUE.")
    }

    lw <- spdep::mat2listw(W, style = style, zero.policy = zero.policy)
    W  <- spdep::listw2mat(lw)
    nb <- lw$neighbours

    diag$kernel_bandwidth <- bandwidth
  }

  # Output
  info <- list(
    method = method,
    contiguity = if (method == "contiguity") contiguity else NULL,
    distance   = if (method == "distance") distance else NULL,
    kernel     = if (method == "kernel") kernel else NULL,
    style      = style,
    lonlat     = lonlat,
    fallback   = if (method == "contiguity") fallback else NULL,
    zero.policy= zero.policy
  )

  if (output == "matrix") return(W)
  if (output == "listw")  return(lw)
  if (output == "nb")     return(nb)

  list(W = W, listw = lw, nb = nb, info = info, diag = diag)
}
