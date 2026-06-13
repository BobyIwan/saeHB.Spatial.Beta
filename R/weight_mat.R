#' @docType data
#' @title Row-Standardized Spatial Weight Matrix
#'
#' @description
#' A row-standardized proximity matrix (\code{W}) generated from a 6x6 regular grid using Queen contiguity.
#' This matrix is mathematically suitable for the Spatial Simultaneous Autoregressive (SAR) model and Moran's I test.
#'
#' @format A \code{36 x 36} numeric matrix. The values are numbers in the interval \code{[0,1]} representing the proximity of the row and column areas. The sum of the values in each row is exactly 1.
#'
#' @usage data(weight_mat)
"weight_mat"
