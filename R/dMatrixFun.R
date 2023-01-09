dMatrixFun = function(indx, n){
  e.vec=matrix(0,n,1)
  e.vec[indx[1],]=1
  e.vec[indx[2],]=(-1)
  return(e.vec)
}
