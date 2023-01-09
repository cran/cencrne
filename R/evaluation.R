#' Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#'
#' @author Mingyang Ren.
#' @usage evaluation(Z.hat, Z.true, cluster.matrix.hat, cluster.matrix.true,
#'                   P.true, Theta.true, K.hat=4, K.true=4)
#'
#' @description The evaluation function for Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#' @param Z.hat A n * r matrix, the estimated embedding vectors corresponding to n nodes.
#' @param Z.true A n * r matrix, the true embedding vectors corresponding to n nodes.
#' @param cluster.matrix.hat A n * n estimated membership matrix, whose (i,j)-element is 1, if nodes i and j are estimated to belong to the same community, and 0, otherwise.
#' @param cluster.matrix.true A n * n true membership matrix, whose (i,j)-element is 1, if nodes i and j belong to the same community, and 0, otherwise.
#' @param P.true A n * n true probability matrix.
#' @param Theta.true A n * n true matrix: Z.true %*% t(Z.true).
#' @param K.hat The true number of communities.
#' @param K.true The estimated number of communities.
#'
#'
#' @return A vector including five evaluation index.
#'         prop. 1: the estimated and actual number of communities are equal; 0: not equal.
#' @export
#'
#'

evaluation = function(Z.hat, Z.true, cluster.matrix.hat, cluster.matrix.true,
                      P.true, Theta.true, K.hat=4, K.true=4){
  if(length(Z.hat) > 0){
    n = dim(Z.hat)[1]
    p = dim(Z.hat)[2]
  } else {
    n = dim(cluster.matrix.true)[1]
  }

  diff.clu.mat = as.matrix(cluster.matrix.hat - cluster.matrix.true)
  group_error_rate = sum(abs(diff.clu.mat[upper.tri(diff.clu.mat)]))/(n*(n-1)/2)

  Theta.hat = Z.hat %*% t(Z.hat)
  P.hat = exp(Theta.hat)/(1+exp(Theta.hat))
  prob.error = sum((P.hat-P.true)^2) / sum((P.true)^2)
  Theta.error = sum((Theta.hat-Theta.true)^2) / sum((Theta.true)^2)

  res = c(as.numeric(K.hat==K.true), K.hat, group_error_rate, prob.error, Theta.error)
  res = as.data.frame(t(res))
  names(res) = c("prop", "K.hat", "cluster.error", "prob.error", "Theta.error")

  return(res)
}

