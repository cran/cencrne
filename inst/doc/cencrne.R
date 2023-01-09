## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(cencrne)
#  # example.data
#  data(example.data)
#  A                   = example.data$A
#  K.true              = example.data$K.true
#  Z.true              = example.data$Z.true
#  B.true              = example.data$B.true
#  P.true              = example.data$P.true
#  Theta.true          = example.data$Theta.true
#  cluster.matrix.true = example.data$cluster.matrix.true
#  
#  n       = dim(A)[1]
#  lam.max = 3
#  lam.min = 0.5
#  lam1.s  = 2/log(n)
#  lam2.s  = sqrt(8*log(n)/n)
#  lam3.s  = 1/8/log(n)/sqrt(n)
#  lambda  = genelambda.obo(nlambda1=3,lambda1_max=lam.max*lam1.s,lambda1_min=lam.min*lam1.s,
#                           nlambda2=10,lambda2_max=lam.max*lam2.s,lambda2_min=lam.min*lam2.s,
#                           nlambda3=1,lambda3_max=lam.max*lam3.s,lambda3_min=lam.min*lam3.s)
#  

## ----eval=FALSE---------------------------------------------------------------
#  sample.index.n = rbind(combn(n,2),1:(n*(n-1)/2))
#  int.list       = gen.int(A)
#  Z.int          = int.list$Z.int
#  B.int          = int.list$B.int
#  res            = network.comm.num(A, sample.index.n, lambda, Z.int, B.int)
#  
#  # output results
#  K.hat = res$Opt_K # the estimated number of communities
#  Z.hat = res$Opt_Z # the estimated embedding vectors corresponding to n nodes
#  cluster.matrix.hat = res$Opt_cluster.matrix # the n * n estimated membership matrix
#  evaluation(Z.hat, Z.true, cluster.matrix.hat, cluster.matrix.true,
#             P.true, Theta.true, K.hat, K.true)
#  

