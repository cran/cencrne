#' Some example data
#'
#' @docType data
#' @name example.data
#' @format A list including:
#' A: An observed n * n adjacency matrix of undirected graph, n=360.
#' K.true: The estimated number of communities.
#' Z.true: A n * r matrix, the true embedding vectors corresponding to n nodes, n=360, r=5.
#' B.true: A n * r matrix, the true community centers corresponding to n nodes.
#' P.true: A n * n true probability matrix.
#' Theta.true: A n * n true matrix: Z.true %*% t(Z.true).
#' cluster.matrix.true: A n * n true membership matrix, whose (i,j)-element is 1, if nodes i and j belong to the same community, and 0, otherwise.
#' @source Simulated data
#' @examples data(example.data)
NULL
