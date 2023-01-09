cluster.memb = function(cluster.matrix){
  n = dim(cluster.matrix)[1]
  memb = rep(0,n)
  num_subgroup = unique(apply(cluster.matrix, 1, function(a){which(a == 1)}))
  if(is.list(num_subgroup)){
    K = length(num_subgroup)
    for (l in 1:K) {
      memb[num_subgroup[[l]]] = l
    }
  } else {
    num_subgroup = unique(t(apply(cluster.matrix, 1, function(a){which(a == 1)})))
    K = dim(num_subgroup)[1]
    for (l in 1:K) {
      memb[num_subgroup[l,]] = l
    }
  }
  return(memb)
}
