#' @title saeHB.Spatial.Beta : Small Area Estimation Hierarchical Bayes for Spatial Beta Model
#'
#' @description Provides several functions and datasets for area-level Small Area Estimation using the Hierarchical Bayesian (HB) method.
#' Model-based estimators are designed for variables of interest that follow a Beta distribution (proportions bounded between 0 and 1).
#' The package supports spatial structures under the Simultaneous Autoregressive (SAR) model and the Leroux Conditional Autoregressive (CAR) model.
#' It also accommodates survey design effect (DEFF) adjustments to handle complex survey data. The \code{rjags} package is employed to obtain parameter estimates via Markov Chain Monte Carlo (MCMC).
#'
#' @section Author(s):
#' Boby Iwan, Cucu Sumarni
#'
#' \strong{Maintainer}: Boby Iwan \email{bobyiwanboby2122@@gmail.com}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{betadeff_sar}}}{Estimates small area means using Spatial SAR Model with Beta distribution and Design Effect (DEFF) adjustments.}
#'   \item{\code{\link{beta_sar}}}{Estimates small area means using Spatial SAR Model with Beta distribution without DEFF adjustments (estimates a global precision parameter).}
#'   \item{\code{\link{betadeff_lerouxcar}}}{Estimates small area means using Spatial Leroux CAR Model with Beta distribution and Design Effect (DEFF) adjustments.}
#'   \item{\code{\link{beta_lerouxcar}}}{Estimates small area means using Spatial Leroux CAR Model with Beta distribution without DEFF adjustments.}
#'   \item{\code{\link{betadeff_nonspatial}}}{Estimates small area means using a Non-Spatial Beta Model with Independent and Identically Distributed (IID) random effects and DEFF adjustments.}
#'   \item{\code{\link{beta_nonspatial}}}{Estimates small area means using a Non-Spatial Beta Model without DEFF adjustments.}
#'   \item{\code{\link{build_w}}}{A utility function to construct spatial weights matrices (contiguity, distance, or kernel) required for spatial modeling.}
#'   \item{\code{\link{moran_test}}}{A diagnostic function to perform Moran's I test for spatial autocorrelation.}
#' }
#'
#' @section Reference:
#' \itemize{
#'   \item{Rao, J. N. K., & Molina, I. (2015). Small Area Estimation (2nd Edition). New Jersey: John Wiley and Sons, Inc. <doi:10.1002/9781118735855>.}
#'   \item{Kubacki, J., & Jedrzejczak, A. (2016). Small Area Estimation of Income Under Spatial SAR Model. Statistics in Transition New Series, Vol. 17, No. 3, pp. 365--390. <doi:10.59170/stattrans-2016-022>.}
#'   \item{Leroux, B. G., Lei, X., & Breslow, N. (2000). Estimation of Disease Rates in Small Areas: A New Mixed Model for Spatial Dependence. In M. E. Halloran & D. Berry (Eds.), Statistical Models in Epidemiology, the Environment, and Clinical Trials (Vol. 116, pp. 179--191). New York: Springer. <doi:10.1007/978-1-4612-1284-3_4>.}
#'   \item{Chung, H. C., & Datta, G. S. (2020). Bayesian Hierarchical Spatial Models for Small Area Estimation. Research Report Series. Washington, D.C.: U.S. Census Bureau.}
#' }
#'
#' @keywords internal
"_PACKAGE"
#'
#' @import rjags
#' @import coda
#' @import stats
#' @import grDevices
#' @import graphics
#' @import sf
#' @import spdep
#'
NULL
