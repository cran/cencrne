network.Z = function(A, Z.int, lambda3=0, a=3, kappa=1, alpha=1,
                     eps = 1e-2, niter = 20, niter.Z=5){
  n = as.integer(dim(A)[1])
  p = dim(Z.int)[2]
  nc = n*(n-1)

  Z.new = Z.int

  iter = 0
  eps.diff.Z = 10
  while( eps.diff.Z >= eps  && iter<niter )
  {
    iter = iter+1
    Z.old = Z.new
    Theta.old = Z.old %*% t(Z.old)
    P.old = exp(Theta.old)/(1+exp(Theta.old))
    P.old[Theta.old > 700] = 1
    diag(P.old) = 0

    # update Z
    alphaZ = 2*alpha
    obj.value.old = 1
    obj.value.new = 10
    iter.Z=0
    while( obj.value.new > obj.value.old  && iter.Z<niter.Z ){
      iter.Z = iter.Z+1
      alphaZ = alphaZ/2
      g.L = - 2/n * t(A - P.old) %*% Z.old
      Z.m = sqrt(apply(Z.old*Z.old, 2, sum))
      g.J3 = t( t(Z.old) * c(mcp_d(Z.m,lambda3) / (Z.m+10^(-10))) )
      Z.new = Z.old - alphaZ*(g.L+g.J3)
      obj.value.old = obj.value.L(A, Z.old)
      obj.value.new = obj.value.L(A, Z.new)
    }

    eps.diff.Z = sqrt(sum((Z.new - Z.old)^2)) / max(sqrt(sum((Z.new)^2)), sqrt(sum((Z.old)^2)) )
  }

  Z.hat = Z.new
  return(list(Z.hat=Z.hat))
}
