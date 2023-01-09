obj.value = function(A, ZZ, BB, lambda1, lambda2, lambda3, Omega.A){
  n = dim(ZZ)[1]
  p = dim(ZZ)[2]

  Theta.hat = ZZ %*% t(ZZ)
  likelihood.mat =  log(1 + exp(Theta.hat)) - Theta.hat * A
  likelihood.mat[Theta.hat > 700] = Theta.hat[Theta.hat > 700] - Theta.hat[Theta.hat > 700] * as.matrix(A)[as.matrix(Theta.hat) > 700]
  L = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1)

  J1 = lambda1*sum((ZZ - BB)^2)

  BB.diff = t(t(BB) %*% Omega.A)
  BB.diff.norm = Matrix::Matrix(sqrt(apply(BB.diff*BB.diff, 1, sum)), sparse = T)
  J2 = sum(MCP_func(BB.diff.norm, lambda2))

  Zm.norm = Matrix::Matrix(sqrt(apply(ZZ*ZZ, 1, sum)), sparse = T)
  J3 = sum(MCP_func(Zm.norm, lambda3))

  obj.val = L + J1 + J2 + J3
  return(obj.val)
}
