#' @docType data
#' @title Row-Standardized Spatial Weight Matrix (SAR)
#'
#' @description
#' A row-standardized proximity matrix (\code{W}) generated from a 6x6 regular grid using Queen contiguity.
#' This matrix is mathematically suitable for the Spatial Simultaneous Autoregressive (SAR) model and Moran's I test.
#'
#' @format A \code{36 x 36} numeric matrix.
#' @usage data(weight_mat)
"weight_mat"
