#' Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#'
#' @author Mingyang Ren.
#' @references Ren, M., Zhang S., Zhang Q. and Ma S. (2022). Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#' @usage gen.int(A, R=8, K.max0=8,rand.seed=123,
#'                lambda3=0, a=3, kappa=1, alpha=1,
#'                eps = 1e-2, niter = 20, niter.Z=5)
#'
#' @description The function generating the initial values.
#' @param A An observed n * n adjacency matrix of undirected graph.
#' @param R Int, the relatively large dimension of embedding vectors given in advance.
#' @param K.max0 The relatively large upper bound of the number of communities given in advance to generate initial values of B.
#' @param rand.seed The random seed of generating initial value.
#' @param lambda3 A float value, the tuning parameter for sparsity of Z.
#' @param a A float value, regularization parameter in MCP, the default setting is 3.
#' @param kappa A float value, the penalty parameter in ADMM algorithm, the default setting is 1.
#' @param alpha A float value, the step size of coordinate descent algorithm updating Z, the default setting is 1.
#' @param eps A float value, algorithm termination threshold.
#' @param niter Int, maximum number of cycles of the overall ADMM algorithm.
#' @param niter.Z Int, maximum number of cycles of coordinate descent algorithm updating Z.
#'
#' @return A list including all estimated parameters and the BIC values with all choices of given tuning parameters, and the selected optional parameters.
#'         Opt_Z: A n * r matrix, the estimated embedding vectors corresponding to n nodes;
#'         Opt_B: A n * r matrix, the estimated community centers corresponding to n nodes;
#'         Opt_K: Int, the estimated number of communities;
#'         Opt_member: A n-dimensional vector, describing the membership of n nodes;
#'         Opt_cluster.matrix: A n * n membership matrix, whose (i,j)-element is 1, if nodes i and j belong to the same community, and 0, otherwise.
#' @export
#'
#'
#' @examples
#' library(cencrne)
#' data(example.data)
#' A                   = example.data$A
#' K.true              = example.data$K.true
#' Z.true              = example.data$Z.true
#' B.true              = example.data$B.true
#' P.true              = example.data$P.true
#' Theta.true          = example.data$Theta.true
#' cluster.matrix.true = example.data$cluster.matrix.true
#'
#' n              = dim(A)[1]
#' sample.index.n = rbind(combn(n,2),1:(n*(n-1)/2))
#' int.list       = gen.int(A)
#' Z.int          = int.list$Z.int
#' B.int          = int.list$B.int
#'
#'
#'

gen.int = function(A, R=8, K.max0=8,rand.seed=123,
                   lambda3=0, a=3, kappa=1, alpha=1,
                   eps = 1e-2, niter = 20, niter.Z=5){
  n = as.integer(dim(A)[1])
  set.seed(rand.seed)
  Z.int = matrix(rnorm(n*R,0,1),n,R)
  Z.int.list = network.Z(A, Z.int, lambda3=lambda3, a=a, kappa=kappa,
                         alpha=alpha, eps = eps, niter = niter, niter.Z=niter.Z)
  Z.int = as.matrix(Z.int.list$Z.hat)
  kmeans.clust = kmeans(Z.int, K.max0)
  memb = kmeans.clust$cluster
  B.int = Z.int
  K.memb = as.data.frame(table(memb))
  for (i in 1:dim(K.memb)[1]) {
    nemb.k = which(memb == as.numeric(K.memb[i,1]))
    if(length(nemb.k) == 1){
      B.int[nemb.k,] = t(matrix(rep(apply(t(Z.int[nemb.k,]),2,mean),length(nemb.k)),ncol=length(nemb.k)))
    } else {
      B.int[nemb.k,] = t(matrix(rep(apply(Z.int[nemb.k,],2,mean),length(nemb.k)),ncol=length(nemb.k)))
    }

  }
  return(list(Z.int=Z.int,B.int=B.int,memb=memb))
}
