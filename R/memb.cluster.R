memb.cluster = function(memb){
  n = length(memb)
  capgfl.matrix = matrix(0,n,n)
  for(i in 1:(n-1)){
    for (j in (i+1):(n)) {
      capgfl.matrix[i,j] <- as.numeric(memb[i] == memb[j])
    }
  }
  capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
  return(capgfl.matrix2)
}
