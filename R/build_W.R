#' @title Build Spatial Weights Matrix
#'
#' @description This function constructs spatial weights matrices (row-standardized \eqn{W} or binary adjacency \eqn{B}) for Hierarchical Bayesian (HB) Beta Spatial modeling under Spatial Autoregressive (SAR) and Leroux Conditional Autoregressive (CAR) models. It supports various methods including Contiguity, Distance-based, and Kernel-based weights, and provides a robust fallback mechanism to automatically connect isolated areas (islands).
#'
#' @details
#' The function supports the following spatial weight construction methods:
#' \itemize{
#'   \item \strong{Contiguity:} Queen, Rook, and Bishop.
#'   \item \strong{Distance-based:} K-Nearest Neighbors (KNN), Inverse Distance, and Exponential.
#'   \item \strong{Kernel-based:} Uniform, Gaussian, Triangular, Epanechnikov, and Quartic.
#' }
#'
#' \strong{Important note on \code{style} for HB Beta Spatial models:}
#' \itemize{
#'   \item \code{style = "W"}: Returns a row-standardized weights matrix (rows sum to 1). This is the required format for HB Beta Spatial Autoregressive (SAR) models.
#'   \item \code{style = "B"}: Returns a binary adjacency matrix (all positive weights are forcefully converted to 1). This is the required format for HB Beta Conditional Autoregressive (CAR) models, specifically the Leroux CAR model.
#' }
#' For continuous weighting methods (Inverse Distance, Exponential, Kernel), using \code{style = "B"} will ignore the calculated continuous weights and binarize the connections. Therefore, if you intend to use Leroux CAR models, Contiguity or KNN are naturally recommended. For SAR models, all methods (including continuous ones) are fully supported and will preserve their weights using \code{style = "W"}.
#'
#' @param data An \code{sf} object containing spatial polygons/points, or a standard data frame.
#' @param coords An \eqn{N \times 2} matrix of spatial coordinates. Required only if \code{data} is not an \code{sf} object.
#' @param method A string indicating the spatial weight construction method. Options are \code{"contiguity"}, \code{"distance"}, or \code{"kernel"}.
#' @param contiguity A string indicating the contiguity type. Options are \code{"queen"}, \code{"rook"}, or \code{"bishop"}. Required if \code{method = "contiguity"}.
#' @param fallback A string indicating the fallback method for isolated areas when using contiguity. Options are \code{"knn"}, \code{"distance"}, or \code{"none"}. Default is \code{"knn"}.
#' @param fallback_k An integer specifying the number of neighbors for the fallback method. Used if \code{fallback = "knn"}. Default is \code{2}.
#' @param fallback_dmax A numeric specifying the maximum distance for the fallback method. Required if \code{fallback = "distance"}.
#' @param distance A string indicating the distance-based type. Options are \code{"knn"}, \code{"inverse_distance"}, or \code{"exponential"}. Required if \code{method = "distance"}.
#' @param k An integer specifying the number of nearest neighbors. Used if \code{distance = "knn"}. Default is \code{2}.
#' @param dmax A numeric specifying the maximum distance threshold. Required if \code{distance} is \code{"inverse_distance"} or \code{"exponential"}.
#' @param power A numeric specifying the decay power. Used if \code{distance = "inverse_distance"}. Default is \code{1}.
#' @param alpha A numeric specifying the decay parameter. Used if \code{distance = "exponential"}. Default is \code{1}.
#' @param epsilon A small numeric value to prevent division by zero in inverse distance calculations. Default is \code{1e-12}.
#' @param kernel A string indicating the type of spatial kernel. Options are \code{"uniform"}, \code{"gaussian"}, \code{"triangular"}, \code{"epanechnikov"}, or \code{"quartic"}. Required if \code{method = "kernel"}.
#' @param bandwidth A numeric specifying the bandwidth (\eqn{h}) for kernel weights. Required if \code{method = "kernel"}.
#' @param lonlat Logical; if \code{TRUE}, coordinates are treated as Longitude/Latitude (degrees), and distances are calculated in kilometers. If \code{FALSE}, distances are calculated using Euclidean geometry in the native units of the coordinates (typically meters for projected CRS). Default is \code{TRUE}.
#' @param style A character string specifying the spatial weights coding scheme. Options are \code{"W"} (row-standardized, for HB Beta SAR models) or \code{"B"} (binary adjacency, for HB Beta Leroux CAR models). Default is \code{"W"}.
#' @param zero.policy Logical; if \code{TRUE}, areas with no neighbors are allowed to have zero-weight rows. Default is \code{TRUE}.
#' @param output A character string specifying the desired format of the returned object. Options are \code{"all"}, \code{"matrix"}, \code{"listw"}, or \code{"nb"}. Default is \code{"all"}.
#'
#' @return Depending on the \code{output} argument, this function returns:
#' \itemize{
#'   \item \strong{\code{"matrix"}}: An \eqn{N \times N} spatial weights matrix (\code{W}). Used as the spatial weight input for SAR and Leroux CAR models.
#'   \item \strong{\code{"listw"}}: A \code{listw} object. Used as the input for \code{moran_test()}.
#'   \item \strong{\code{"nb"}}: An \code{nb} (neighborhood) object representing the list of neighbors for each area.
#'   \item \strong{\code{"all"}}: A comprehensive list containing \code{W}, \code{listw}, \code{nb}, \code{info}, and \code{diag}.
#' }
#'
#' @examples
#' library(sf)
#'
#' # 1. Contiguity Method (Requires sf polygons)
#' # Create a simple 3x3 polygon grid
#' bbox <- st_bbox(c(xmin = 0, ymin = 0, xmax = 3, ymax = 3))
#' grid <- st_make_grid(bbox, n = c(3, 3))
#' grid_sf <- st_sf(id = 1:9, geometry = grid)
#'
#' # Build spatial weights with output = "all"
#' W_obj <- build_w(
#'   data = grid_sf,
#'   method = "contiguity",
#'   contiguity = "queen",
#'   style = "W",
#'   output = "all"
#' )
#'
#' # Extract outputs from the list
#' W_mat   <- W_obj$W      # Spatial weights matrix (for SAR/CAR)
#' lw_obj  <- W_obj$listw  # listw object (for moran_test)
#' W_diag  <- W_obj$diag   # Diagnostic info (isolated areas, etc.)
#' head(W_mat)
#'
#' # Setup Coordinates for Distance & Kernel
#' # Generate random Longitude and Latitude coordinates for 10 areas
#' set.seed(123)
#' lon <- runif(10, min = 100, max = 140)
#' lat <- runif(10, min = -10, max = 10)
#' coords <- cbind(lon, lat)
#'
#' # 2. Distance Method (Using coordinates)
#' # Build row-standardized KNN weights (style = "W") for SAR models
#' W_knn <- build_w(
#'   data = NULL,
#'   coords = coords,
#'   method = "distance",
#'   distance = "knn",
#'   k = 2,
#'   lonlat = TRUE,
#'   style = "W",
#'   output = "matrix"
#' )
#' head(W_knn)
#'
#' # 3. Kernel Method (Using coordinates)
#' # Build binary adjacency Kernel weights (style = "B") for Leroux CAR models
#' W_kernel <- build_w(
#'   data = NULL,
#'   coords = coords,
#'   method = "kernel",
#'   kernel = "gaussian",
#'   bandwidth = 500,
#'   lonlat = TRUE,
#'   style = "B",
#'   output = "matrix"
#' )
#' head(W_kernel)
#'
#' @import sf
#' @import spdep
#' @importFrom stats na.omit
#'
#' @export build_w
build_w <- function(
    data,
    coords = NULL,
    method = c("contiguity", "distance", "kernel"),
    contiguity = c("queen", "rook", "bishop"),
    fallback = c("knn", "distance", "none"),
    fallback_k = 2,
    fallback_dmax = NULL,
    distance = c("knn", "inverse_distance", "exponential"),
    k = 2,
    dmax = NULL,
    power = 1,
    alpha = 1,
    epsilon = 1e-12,
    kernel = c("uniform", "gaussian", "triangular", "epanechnikov", "quartic"),
    bandwidth = NULL,
    lonlat = TRUE,
    style = c("W", "B"),
    zero.policy = TRUE,
    output = c("all", "matrix", "listw", "nb")
) {

  method     <- match.arg(method)
  contiguity <- match.arg(contiguity)
  fallback   <- match.arg(fallback)
  distance   <- match.arg(distance)
  kernel     <- match.arg(kernel)
  style      <- match.arg(style)
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

      diag$k <- k
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
          stop("Some areas have no neighbors with current dmax. Increase dmax.")
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

      if (style == "B") {
        lw <- spdep::nb2listw(nb, style = "B", zero.policy = zero.policy)
      } else {
        lw <- spdep::nb2listw(nb, glist = glist, style = "W", zero.policy = zero.policy)
      }

      W <- spdep::listw2mat(lw)
      diag$used_dmax <- !is.null(dmax)
      if (!is.null(dmax)) {
        diag$dmax <- dmax
      } else {
        diag$k <- k
      }
      diag$isolates_after <- which(spdep::card(nb) == 0)
      diag$n_isolates_after <- length(diag$isolates_after)
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
      diag(W) <- 0
    }

    if (any(rowSums(W) == 0) && !zero.policy) {
      stop("Some areas have no neighbors with current bandwidth. Increase bandwidth (h).")
    }

    if (style == "B") {
      W[W > 0] <- 1
      lw <- spdep::mat2listw(W, style = "B", zero.policy = zero.policy)
    } else {
      lw <- spdep::mat2listw(W, style = "W", zero.policy = zero.policy)
    }

    W  <- spdep::listw2mat(lw)
    nb <- lw$neighbours

    diag$kernel_bandwidth <- bandwidth
    diag$isolates_after <- which(rowSums(W) == 0)
    diag$n_isolates_after <- length(diag$isolates_after)
  }

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
