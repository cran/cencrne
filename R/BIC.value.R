BIC.value = function(A, Z.hat, K.hat, adjust=FALSE, B.hat, ad.B=FALSE){
  n = dim(Z.hat)[1]
  p = dim(Z.hat)[2]
  r.hat = 1

  Theta.hat = Z.hat %*% t(Z.hat)
  likelihood.mat =  log(1 + exp(Theta.hat)) - Theta.hat * A
  if(ad.B){
    fitness = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1) + sum((Z.hat-B.hat)^2)/n
  } else {
    fitness = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1)
  }

  if(adjust){Cn = log(n*p)} else {Cn = 1}
  degree = Cn*r.hat*K.hat*log(n)  / n
  BICvalue = fitness + degree
  return(list(BIC=BICvalue, fitness=fitness, degree=degree))
}
