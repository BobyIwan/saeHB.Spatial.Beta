#' @title saeHB.Spatial.Beta : Small Area Estimation Hierarchical Bayes for Spatial Beta Model
#'
#' @description Provides several functions and datasets for area-level Small Area Estimation using the Hierarchical Bayesian (HB) method.
#' Model-based estimators are designed for variables of interest that follow a Beta distribution (proportions bounded between 0 and 1).
#' The package supports both non-spatial models and spatial models based on the Simultaneous Autoregressive (SAR) and Leroux Conditional Autoregressive (CAR) structures, with optional survey design effect (DEFF) adjustments.
#' In addition, it provides utility functions for constructing spatial weights matrices and performing spatial autocorrelation diagnostics.
#' The \code{rjags} package is used to obtain posterior estimates via Markov Chain Monte Carlo (MCMC).
#'
#' @section Author(s):
#' Boby Iwan, Cucu Sumarni
#'
#' \strong{Maintainer}: Boby Iwan \email{bobyiwanboby2122@@gmail.com}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{betadeff_sar}}}{Estimates small area proportions using a Hierarchical Bayes (HB) method under a Spatial SAR Model with a Beta distribution, incorporating survey design effect (DEFF) adjustments.}
#'   \item{\code{\link{beta_sar}}}{Estimates small area proportions using a Hierarchical Bayes (HB) method under a Spatial SAR Model with a Beta distribution without DEFF adjustments, by estimating the unknown precision parameter.}
#'   \item{\code{\link{betadeff_lerouxcar}}}{Estimates small area proportions using a Hierarchical Bayes (HB) method under a Spatial Leroux CAR Model with a Beta distribution, incorporating survey design effect (DEFF) adjustments.}
#'   \item{\code{\link{beta_lerouxcar}}}{Estimates small area proportions using a Hierarchical Bayes (HB) method under a Spatial Leroux CAR Model with a Beta distribution without DEFF adjustments, by estimating the unknown precision parameter.}
#'   \item{\code{\link{betadeff_nonspatial}}}{Estimates small area proportions using a Hierarchical Bayes (HB) method under a Non-Spatial Model with a Beta distribution and Independent and Identically Distributed (IID) random effects, incorporating DEFF adjustments.}
#'   \item{\code{\link{beta_nonspatial}}}{Estimates small area proportions using a Hierarchical Bayes (HB) method under a Non-Spatial Model with a Beta distribution and IID random effects without DEFF adjustments, by estimating the unknown precision parameter.}
#'   \item{\code{\link{build_w}}}{A utility function to construct spatial weights matrices (contiguity, distance, or kernel) required for spatial modeling.}
#'   \item{\code{\link{moran_test}}}{A diagnostic function to perform Moran's I test for spatial autocorrelation.}
#' }
#'
#' @section Reference:
#' \itemize{
#'   \item{Rao, J. N. K., & Molina, I. (2015). Small Area Estimation (2nd Edition). New Jersey: John Wiley and Sons, Inc. <doi:10.1002/9781118735855>.}
#'   \item{Liu, B. (2009). Hierarchical Bayes estimation and empirical best prediction of small-area proportions. <https://api.drum.lib.umd.edu/server/api/core/bitstreams/cb8e2cbf-441e-4f0f-b4b3-6182f3cf24de/content>.}
#'   \item{Liu, B., Lahiri, P., & Kalton, G. (2014). Hierarchical Bayes Modeling of Survey-Weighted Small Area Proportions. Statistics Canada. <https://www150.statcan.gc.ca/n1/pub/12-001-x/2014001/article/14030-eng.pdf>.}
#'   \item{Kubacki, J., & Jedrzejczak, A. (2016). Small Area Estimation of Income Under Spatial SAR Model. Statistics in Transition New Series, Vol. 17, No. 3, pp. 365--390. <doi:10.59170/stattrans-2016-022>.}
#'   \item{Leroux, B. G., Lei, X., & Breslow, N. (2000). Estimation of Disease Rates in Small Areas: A New Mixed Model for Spatial Dependence. In M. E. Halloran & D. Berry (Eds.), Statistical Models in Epidemiology, the Environment, and Clinical Trials (Vol. 116, pp. 179--191). New York: Springer. <doi:10.1007/978-1-4612-1284-3_4>.}
#'   \item{Chung, H. C., & Datta, G. S. (2020). Bayesian Hierarchical Spatial Models for Small Area Estimation. Research Report Series. Washington, D.C.: U.S. Census Bureau. <https://www.census.gov/content/dam/Census/library/working-papers/2020/adrm/RRS2020-07.pdf>.}
#'   \item{Anselin, L. (1988). Spatial Econometrics: Methods and Models. Dordrecht: Springer Netherlands. <doi:10.1007/978-94-015-7799-1>.}
#'   \item{Anselin, L., & Morrison, S. (2019). Spatial Weights as Distance Functions. <https://spatialanalysis.github.io/lab_tutorials/Spatial_Weights_as_Distance_Functions.html>.}
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
