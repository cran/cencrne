#' Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#'
#' @author Mingyang Ren, Sanguo Zhang, Junhui Wang. Maintainer: Mingyang Ren <renmingyang17@mails.ucas.ac.cn>.
#' @references Ren, M., Zhang S. and Wang J. (2022). Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#' @usage network.comm.num(A, sample.index.n, lambda, Z.int, B.int,
#'                         a=3, kappa=1, alpha=1, eps=5e-2, niter=20,
#'                niter.Z=5, update.B="ADMM",local.oppro=FALSE, merge.all=TRUE,
#'                ad.BIC=FALSE, Fully.Connected=TRUE, trace=FALSE,
#'                line.search=TRUE, ad.BIC.B=FALSE)
#'
#' @description The main function for Consistent Estimation of the Number of Communities via Regularized Network Embedding.
#' @param A An observed n * n adjacency matrix of undirected graph.
#' @param sample.index.n A 3 * (n*(n-1)/2) matrix, all pairs of integers from 1 to n.
#' @param lambda A list, the sequences of the tuning parameters (\eqn{\lambda_1}{\lambda_1}, \eqn{\lambda_2}{\lambda_2}, and \eqn{\lambda_3}{\lambda_3}).
#' @param Z.int A n * r matrix, the initial values of embedding vectors corresponding to n nodes.
#' @param B.int A n * r matrix, the initial values of community centers corresponding to n nodes.
#' @param a A float value, regularization parameter in MCP, the default setting is 3.
#' @param kappa A float value, the penalty parameter in ADMM algorithm, the default setting is 1.
#' @param alpha A float value, the step size of coordinate descent algorithm updating Z, the default setting is 1.
#' @param eps A float value, algorithm termination threshold.
#' @param niter Int, maximum number of cycles of the overall ADMM algorithm.
#' @param niter.Z Int, maximum number of cycles of coordinate descent algorithm updating Z.
#' @param update.B The optimization algorithm updating B, which can be selected "ADMM" (default) and "AMA".
#' @param local.oppro The logical variable, whether to use local approximations when updating Z, the default setting is F.
#' @param Fully.Connected Whether to use the all pairs (i,j) in fusion penalty, the default setting is T. If F, the pairs (i,j) in fusion penalty will be determined by the observed n * n adjacency matrix A.
#' @param trace Whether to output the intermediate process of the algorithm.
#' @param line.search Linear search or not, the default setting is T.
#' @param ad.BIC Whether to use the adjusted BIC, the default setting is F.
#' @param merge.all Whether to merge pairs of nodes indirectly connected (but without the direct edge) in the estimated community membership matrix.
#' @param ad.BIC.B Whether the BIC criterion contains terms involving the B matrix, the default setting is F.
#'
#' @return A list including all estimated parameters and the BIC values with all choices of given tuning parameters, and the selected optional parameters.
#'         Opt_Z: A n * r matrix, the estimated embedding vectors corresponding to n nodes;
#'         Opt_B: A n * r matrix, the estimated community centers corresponding to n nodes;
#'         Opt_K: Int, the estimated number of communities;
#'         Opt_member: A n-dimensional vector, describing the membership of n nodes;
#'         Opt_cluster.matrix: A n * n membership matrix, whose (i,j)-element is 1, if nodes i and j belong to the same community, and 0, otherwise.
#'
#'
#' @export
#' @importFrom stats kmeans median rnorm
#' @import Matrix
#'
#' @examples
#' \donttest{
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
#' n       = dim(A)[1]
#' lam.max = 3
#' lam.min = 0.5
#' lam1.s  = 2/log(n)
#' lam2.s  = sqrt(8*log(n)/n)
#' lam3.s  = 1/8/log(n)/sqrt(n)
#' lambda  = genelambda.obo(nlambda1=3,lambda1_max=lam.max*lam1.s,lambda1_min=lam.min*lam1.s,
#'                          nlambda2=10,lambda2_max=lam.max*lam2.s,lambda2_min=lam.min*lam2.s,
#'                          nlambda3=1,lambda3_max=lam.max*lam3.s,lambda3_min=lam.min*lam3.s)
#'
#' sample.index.n = rbind(combn(n,2),1:(n*(n-1)/2))
#' int.list       = gen.int(A)
#' Z.int          = int.list$Z.int
#' B.int          = int.list$B.int
#' res            = network.comm.num(A, sample.index.n, lambda, Z.int, B.int)
#'
#' K.hat = res$Opt_K # the estimated number of communities
#' Z.hat = res$Opt_Z # the estimated embedding vectors corresponding to n nodes
#' cluster.matrix.hat = res$Opt_cluster.matrix # the n * n estimated membership matrix
#' evaluation(Z.hat, Z.true, cluster.matrix.hat, cluster.matrix.true,
#'            P.true, Theta.true, K.hat, K.true)
#' }
#'

network.comm.num = function(A, sample.index.n, lambda, Z.int, B.int,
                            a=3, kappa=1, alpha=1, eps=5e-2, niter=20,
                            niter.Z=5, update.B="ADMM",local.oppro=FALSE, merge.all=TRUE,
                            ad.BIC=FALSE, Fully.Connected=TRUE, trace=FALSE,
                            line.search=TRUE, ad.BIC.B=FALSE){
  n  = as.integer(dim(A)[1])
  p  = dim(Z.int)[2]

  Omega.list     = gen.weight(A, sample.index.n, Fully.Connected=Fully.Connected)
  Omega.A        = Omega.list$Omega.A
  sample.index.A = Omega.list$sample.index.A
  gamma.index    = Omega.list$gamma.index

  if(line.search){
    lam1     = lambda$lambda1
    lam2     = lambda$lambda2
    lam3     = lambda$lambda3
    L1       = length(lam1)
    L2       = length(lam2)
    L3       = length(lam3)
    L        = L1+L2+L3
    aBIC     = rep(0,L)
    K.list   = rep(0,L)
    Z.list   = list()
    B.list   = list()
    V.list   = list()
    BIC.list = list()
    lam.list = list()
    member.list = list()
    cluster.matrix.list = list()
    if(L == 3){
      l=1
      aBIC = rep(0,l)
      K.list  = rep(0,l)
      lambda1 = lam1; lambda2 = lam2; lambda3 = lam3
      lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
      res.list = fusion.network(A, Omega.A, sample.index.A, sample.index.n, gamma.index,
                                lambda1, lambda2, lambda3, Z.int, B.int,
                                a=a, kappa=kappa, alpha=alpha, eps = eps,
                                niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                update.B=update.B, local.oppro=local.oppro,
                                merge.all=merge.all, ad.BIC.B=ad.BIC.B)
      Z.list[[l]]   = res.list$Z.hat
      B.list[[l]]   = res.list$B.hat
      V.list[[l]]   = res.list$V.hat
      BIC.list[[l]] = res.list$BIC.vec
      aBIC[l]       = res.list$BIC.vec$BIC
      K.list[l]     = res.list$K.hat
      lam.list[[l]] = lam.all
      member.list[[l]] = res.list$cluster.results$member
      cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
      if(trace){
        message(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                paste("K =",as.numeric(res.list$K.hat))))
      }
    } else {
      # search lambda2
      lambda1 = 1;lambda3 = median(lam3)
      for (l in 1:L2) {
        lambda2 = lam2[l]

        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, sample.index.n, gamma.index,
                                  lambda1, lambda2, lambda3, Z.int, B.int,
                                  a=a, kappa=kappa, alpha=alpha, eps = eps,
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro,
                                  merge.all=merge.all, ad.BIC.B=ad.BIC.B)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          message(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
      }
      aBIC[is.na(aBIC)] = 10^10
      n_lam2 = which(aBIC[1:L2] == min(aBIC[1:L2]))[1];lambda2 = lam2[n_lam2]

      # search lam1
      for (l1 in 1:L1) {
        lambda1 = lam1[l1];l = L2+l1

        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, sample.index.n, gamma.index,
                                  lambda1, lambda2, lambda3, Z.int, B.int,
                                  a=a, kappa=kappa, alpha=alpha, eps = eps,
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro,
                                  merge.all=merge.all, ad.BIC.B=ad.BIC.B)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          message(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
      }
      aBIC[is.na(aBIC)] = 10^10
      n_lam1 = which(aBIC[(L2+1):(L2+L1)] == min(aBIC[(L2+1):(L2+L1)]))[1];lambda1 = lam1[n_lam1]

      # search lam3
      for (l3 in 1:L3) {
        lambda3 = lam3[l3];l = L2+L1+l3

        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, sample.index.n, gamma.index,
                                  lambda1, lambda2, lambda3, Z.int, B.int,
                                  a=a, kappa=kappa, alpha=alpha, eps = eps,
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro,
                                  merge.all=merge.all, ad.BIC.B=ad.BIC.B)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          message(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
      }
      aBIC[is.na(aBIC)] = 10^10
      n_lam3 = which(aBIC[(L2+L1+1):(L2+L1+L3)] == min(aBIC[(L2+L1+1):(L2+L1+L3)]))[1];lambda3 = lam3[n_lam3]

    }
  } else {
    lam1     = lambda$lambda1
    lam2     = lambda$lambda2
    lam3     = 0
    L1       = length(lam1)
    L2       = length(lam2)
    L3       = length(lam3)
    L        = L1*L2
    aBIC     = rep(0,L)
    K.list   = rep(0,L)
    Z.list   = list()
    B.list   = list()
    V.list   = list()
    BIC.list = list()
    lam.list = list()
    member.list = list()
    cluster.matrix.list = list()

    l = 1
    for (l2 in 1:L2) {
      for (l1 in 1:L1) {
        lambda1 = lam1[l1]
        lambda2 = lam2[l2]
        lambda3 = 0

        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, sample.index.n, gamma.index,
                                  lambda1, lambda2, lambda3, Z.int, B.int,
                                  a=a, kappa=kappa, alpha=alpha, eps = eps,
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro,
                                  merge.all=merge.all, ad.BIC.B=ad.BIC.B)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          message(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
        l = l+1
      }
    }
    aBIC[is.na(aBIC)] = 10^10
  }


  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)
  Opt_lambda = lam.list[[n_lam]]
  Opt_K = K.list[n_lam]
  Opt_Z = Z.list[[n_lam]]
  Opt_B = B.list[[n_lam]]
  Opt_member = member.list[[n_lam]]
  Opt_cluster.matrix = cluster.matrix.list[[n_lam]]
  result = list(BIC=aBIC, K.vec=K.list, Z.list=Z.list, B.list=B.list,
                BIC.list=BIC.list, lam.list=lam.list,
                member.list=member.list, cluster.matrix.list=cluster.matrix.list,
                Opt_num=n_lam, Opt_Z=Opt_Z, Opt_B=Opt_B,
                Opt_K=Opt_K, Opt_member=Opt_member,
                Opt_aBIC=Opt_aBIC, Opt_lambda=Opt_lambda,
                Opt_cluster.matrix=Opt_cluster.matrix)
  return(result)

}


