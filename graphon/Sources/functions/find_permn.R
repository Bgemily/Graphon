
### Find the optimal permutation
find_permn = function(center_cdf_array_from, center_cdf_array_to){
  if (!identical(dim(center_cdf_array_from), dim(center_cdf_array_to))) {
    stop("dim(center_cdf_array_from) and dim(center_cdf_array_to) should be the same.")
  }
  
  N_clus = dim(center_cdf_array_from)[1]
  permn_list = combinat::permn(1:N_clus)
  
  min_dist = Inf
  for (permn in permn_list) {
    
    ### V1: l2 distance
    # dist = sum((center_cdf_array_from[permn, permn, ] - center_cdf_array_to)^2)
    ### V2: cross-correlation
    dist = 0
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        center_cdf_1 = center_cdf_array_from[permn, permn, ][q,k,]
        center_cdf_2 = center_cdf_array_to[q,k,]
        if (var(center_cdf_1)==0) {
          center_cdf_1 = jitter(center_cdf_1)
        }
        if (var(center_cdf_2)==0) {
          center_cdf_2 = jitter(center_cdf_2)
        }
        tmp = 1-max(ccf(x = center_cdf_1, y = center_cdf_2, plot=FALSE)$acf)
        dist = dist + tmp
      }
    }
  
    ########
    
    if(dist < min_dist){
      the_permn = permn
      min_dist = dist
    }
    
  }
  
  return(list(permn = the_permn, dist = min_dist))
}
