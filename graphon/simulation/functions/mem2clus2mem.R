

# transfer vector of membership to list of clusters 
mem2clus = function(membership){
  N_clus = length(unique(membership))
  clusters = vector("list", N_clus)
  for (l in 1:N_clus) {
    clusters[[l]] = which(membership==unique(membership)[l])
  }
  return(clusters)
}


clus2mem = function(clusters){
  membership = unlist(clusters)
  N_clus = length(clusters)
  for (l in 1:N_clus) {
    membership[clusters[[l]]] = l
  }
  return(membership)
}
