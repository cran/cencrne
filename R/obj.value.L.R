obj.value.L = function(A, ZZ){
  n = dim(ZZ)[1]
  p = dim(ZZ)[2]

  Theta.hat = ZZ %*% t(ZZ)
  likelihood.mat =  log(1 + exp(Theta.hat)) - Theta.hat * A
  likelihood.mat[Theta.hat > 700] = Theta.hat[Theta.hat > 700] - Theta.hat[Theta.hat > 700] * as.matrix(A)[as.matrix(Theta.hat) > 700]
  L = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1)

  obj.val = L
  return(obj.val)
}
