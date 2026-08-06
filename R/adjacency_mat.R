#' @docType data
#' @title Binary Adjacency Matrix
#'
#' @description
#' A binary adjacency matrix (\code{B}) generated from a 6x6 regular grid using Queen contiguity.
#' This matrix is mathematically suitable for the Hierarchical Bayesian (HB) Spatial Beta-Leroux CAR model.
#'
#' @format A \code{36 x 36} numeric matrix. The elements take a value of \code{1} if two areas share a common border (are neighbors), and \code{0} otherwise. All diagonal elements are \code{0}.
#'
#' @usage data(adjacency_mat)
"adjacency_mat"
