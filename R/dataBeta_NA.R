#' @docType data
#' @title Synthetic Data with Missing Values for Small Area Estimation
#'
#' @description
#' A synthetic dataset identical to \code{\link{dataBeta}}, but contains 5 missing values (NA) in the
#' variable of interest (\code{y}) to demonstrate the prediction capability of the models for non-sampled areas.
#'
#' @format A data frame with 100 rows and 6 columns:
#' \describe{
#'   \item{domain}{Area ID/name}
#'   \item{y}{Direct estimates of the proportion/variable of interest (0 < y < 1). Contains \code{NA} values.}
#'   \item{x1}{Auxiliary variable 1}
#'   \item{x2}{Auxiliary variable 2}
#'   \item{n_i}{Sample size for each area}
#'   \item{deff}{Survey design effect for each area}
#' }
#'
#' @usage data(dataBeta_NA)
"dataBeta_NA"
