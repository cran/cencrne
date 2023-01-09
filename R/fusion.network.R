fusion.network = function(A, Omega.A, sample.index.A, sample.index.n, gamma.index,
                          lambda1, lambda2, lambda3, Z.int, B.int, a=3, kappa=1, alpha=1,
                          eps = 5e-2, niter = 20, niter.Z=5, update.B="ADMM",
                          local.oppro=F, ad.BIC=F, merge.all=T, ad.BIC.B=F){
  n = as.integer(dim(A)[1])
  p = dim(Z.int)[2]
  nc = n*(n-1)
  K.int = length(unique(apply(B.int*B.int,1,sum)))

  gamma.int = Matrix::Matrix(0,nrow=dim(Omega.A)[2],ncol=p,sparse = T)
  V.int = t(t(B.int) %*% Omega.A)
  Z.new = Z.int
  B.new = B.int
  gamma.new = gamma.int
  V.new = V.int

  iter = 0
  K.hat = K.int
  niter.K = 3
  iter.K = 0
  while( K.hat == K.int && iter.K<niter.K ){
    iter.K = iter.K+1
    eps.diff.Z = 10
    eps.diff.V = 10

    while( eps.diff.Z+eps.diff.V >= eps && iter<niter )
    {

      iter = iter+1
      Z.old = Z.new
      B.old = B.new
      gamma.old = gamma.new
      V.old = V.new

      Theta.old = Z.old %*% t(Z.old)
      P.old = exp(Theta.old)/(1+exp(Theta.old))
      P.old[Theta.old > 700] = 1
      diag(P.old) = 0

      # update Z
      if(!local.oppro){
        alphaZ = 2*alpha
        obj.value.old = 1
        obj.value.new = 10
        iter.Z=0
        while( obj.value.new > obj.value.old  && iter.Z<niter.Z ){
          iter.Z = iter.Z+1
          alphaZ = alphaZ/2
          g.L = - 2/n * t(A - P.old) %*% Z.old
          g.J1 = 2*lambda1 * (Z.old - B.old)
          Z.m = sqrt(apply(Z.old*Z.old, 2, sum))
          g.J3 = t( t(Z.old) * c(mcp_d(Z.m,lambda3) / (Z.m+10^(-10))) )
          Z.new = Z.old - alphaZ*(g.L+g.J1+g.J3)
          obj.value.old = obj.value(A, Z.old, B.old, lambda1, lambda2, lambda3, Omega.A)
          obj.value.new = obj.value(A, Z.new, B.new, lambda1, lambda2, lambda3, Omega.A)
        }
      } else {
        Z.m = apply(Z.old, 2, function(x) sqrt(sum(x^2)))
        c.Z1 = (lambda1 + mcp_d(Z.m,lambda3) / (Z.m+10^(-10)) / 2)^(-1)
        Z.new = t(t( lambda1 * B.old + 1/nc * t(A - P.old) %*% Z.old ) * c.Z1)
      }

      obj.value.old.B = obj.value(A, Z.new, B.old, lambda1, lambda2, lambda3, Omega.A)
      # update B
      if(update.B=="AMA"){
        gamma.up = Matrix::Matrix(t(t(gamma.old) %*% gamma.index),sparse = T)
        B.new = Z.new + gamma.up
      } else if(update.B=="ADMM"){
        kappa.lam = kappa/2/(lambda1+10^(-10))
        gamma.up = Matrix::Matrix(t(t(gamma.old + kappa.lam*V.old) %*% gamma.index),sparse = T)
        Z.new.mean = n*kappa.lam*apply(Z.new, 2, mean)/(1+n*kappa.lam)
        Y.new = 1/(1+n*kappa.lam)*(Z.new + gamma.up)
        B.new = t(t(Y.new) + Z.new.mean)
      } else {warning("The method of updating B must be AMA or ADMM.")}

      # update V
      V.old = t(t(B.new) %*% Omega.A)
      delta = V.old - gamma.old/kappa
      V.new = Matrix::Matrix(t(apply(delta, 1, MCP_soft, lambda=lambda2, kappa=kappa)),sparse = T)
      # update gamma
      gamma.new = gamma.old + kappa*(V.new - V.old)

      eps.diff.Z = sqrt(sum((Z.new - Z.old)^2)/length(Z.old))
      eps.diff.V = sqrt(sum((V.new - V.old)^2)/length(V.old))

      obj.value.new.B = obj.value(A, Z.new, B.new, lambda1, lambda2, lambda3, Omega.A)

    }

    Z.hat = Z.new
    B.hat = B.new
    V.hat = V.new

    V.norm = apply(V.hat*V.hat,1,sum)
    diff.Bij = rep(1,nc/2)
    diff.Bij[sample.index.A[3,V.norm == 0]] = 0
    cluster.results = get.cluster(diff.Bij, n, sample.index.n, merge.all=merge.all)
    K.hat = cluster.results$K.hat
    BIC.vec = BIC.value(A, Z.hat, K.hat, adjust=ad.BIC, B.hat, ad.B=ad.BIC.B)
  }

  Z.hat = Z.new
  B.hat = B.new
  V.hat = V.new

  V.norm = apply(V.hat*V.hat,1,sum)
  diff.Bij = rep(1,nc/2)
  diff.Bij[sample.index.A[3,V.norm == 0]] = 0
  cluster.results = get.cluster(diff.Bij, n, sample.index.n, merge.all=merge.all)
  K.hat = cluster.results$K.hat
  BIC.vec = BIC.value(A, Z.hat, K.hat, adjust=ad.BIC, B.hat)

  return(list(Z.hat=Z.hat, B.hat=B.hat, V.hat=V.hat, K.hat=K.hat,
              cluster.results=cluster.results, BIC.vec=BIC.vec))
}
