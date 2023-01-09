gen.weight = function(A, sample.index.n, Fully.Connected=F){
  n = dim(A)[1]
  A.con = A
  if(Fully.Connected){A.con[A==0] = 1}
  omega.n = NULL
  for (i in 1:(n-1)) {
    index.i = as.matrix(sample.index.n[,sample.index.n[1,] == i])
    omega.n = c(omega.n,index.i[3,match(c(i+which(A.con[i,(i+1):n] == 1)),index.i[2,])])
  }
  Omega = Matrix::Matrix(apply(sample.index.n[1:2,],2,dMatrixFun,n=n),sparse = T)
  Omega.A = Omega[,omega.n]
  sample.index.A = rbind(sample.index.n[,omega.n],1:length(omega.n))

  gamma.index = Matrix::Matrix(0,length(omega.n),n, sparse = T)
  for (i in 1:n) {
    nn.ga.after = which(sample.index.A[1,] == i)
    nn.ga.before = which(sample.index.A[2,] == i)
    if(length(nn.ga.before)>0 & length(nn.ga.after)>0){
      gamma.index[nn.ga.after,i] = 1
      gamma.index[nn.ga.before,i] = -1
    }
    if(i==1){gamma.index[nn.ga.after,i] = 1}
    if(i==n){gamma.index[nn.ga.before,i] = -1}
  }

  return(list(Omega.A=Omega.A, sample.index.A=sample.index.A, gamma.index=gamma.index))
}

