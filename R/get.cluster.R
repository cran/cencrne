get.cluster = function(diff.Bij, n, sample.index.n, merge.all=T){
  size = n
  diff.gfl=diff.Bij
  capgfl.matrix=Matrix::Matrix(0,nrow=size,ncol=size,sparse = T)

  if(length(which(diff.gfl==0))==0){
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = size
    memb = 1:n
  }else if(length(which(diff.gfl==0))==size*(size-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)]=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = 1
    memb = rep(1,n)
  }else{
    sample.index.gfl=sample.index.n[1:2,which(diff.gfl==0)]
    if(length(which(diff.gfl==0))==1){
      capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]]=1
    }else{
      for(i in 1:length(which(diff.gfl==0))){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]]=1
      }
    }
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(as.matrix(capgfl.matrix2)))


    cap = capgfl.matrix2
    cap.cap = apply(cap, 1, function(a){which(a == 1)})
    if(is.list(cap.cap)){num_subgroup = unique(cap.cap)}
    if(is.matrix(cap.cap)){
      uni.cap = unique(t(cap.cap))
      num_subgroup = list()
      for (i in 1:dim(uni.cap)[1]) {
        num_subgroup[[i]] = as.numeric(uni.cap[i,])
      }
    }
    non_inter_list = list()
    vv = 1
    non_inter = c(1:length(num_subgroup))
    repeat{
      a = num_subgroup[[non_inter[1]]]
      KK_k = setdiff(non_inter,non_inter[1])
      non_inter = c()
      i=1
      for (k2 in KK_k) {
        if(length(intersect(a,num_subgroup[[k2]])) > 0){
          a = union(a,num_subgroup[[k2]])
        } else {
          non_inter[i] = k2
          i=i+1
        }
      }
      non_inter_list[[vv]] = a
      vv = vv+1
      if(length(non_inter) == 0){break}
    }

    if(merge.all){
      for (i in 1:dim(cap)[1]) {
        for (k in 1:length(non_inter_list)) {
          if(length(intersect(which(cap[i,]==1),non_inter_list[[k]])) > 0){
            cap[i,non_inter_list[[k]]] = 1
          }
        }
      }
      capgfl.matrix2 = cap
      group.num.gf = nrow(unique(as.matrix(capgfl.matrix2)))
    }

    memb = rep(0,n)
    for (k in 1:length(non_inter_list)) {
      memb[non_inter_list[[k]]] = k
    }

  }

  return(list(K.hat=group.num.gf, cluster.matrix=capgfl.matrix2, member=memb))
}
