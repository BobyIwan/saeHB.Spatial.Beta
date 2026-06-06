#' @docType data
#' @title Synthetic Data for Small Area Estimation using Spatial Beta Model
#'
#' @description
#' A synthetic dataset generated for testing and tutorial purposes of the \code{saeHB.Spatial.Beta} package.
#' The data is generated under a Spatial Simultaneous Autoregressive (SAR) process with a Beta distribution,
#' accommodating survey design effects (DEFF).
#'
#' @format A data frame with 100 rows and 6 columns:
#' \describe{
#'   \item{domain}{Area ID/name}
#'   \item{y}{Direct estimates of the proportion/variable of interest (0 < y < 1)}
#'   \item{x1}{Auxiliary variable 1 (Normal distribution)}
#'   \item{x2}{Auxiliary variable 2 (Normal distribution)}
#'   \item{n_i}{Sample size for each area}
#'   \item{deff}{Survey design effect for each area}
#' }
#'
#' @usage data(dataBeta)
"dataBeta"
