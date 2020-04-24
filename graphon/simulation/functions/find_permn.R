

# find the optimal permutation
find_permn = function(pdf_array_1, pdf_array_2, t_unit = 0.05){
  if (dim(pdf_array_1)[1]!=dim(pdf_array_2)[1]) {
    stop("Shapes of pdf_array_1 and pdf_array_2 should be the same!")
  }
  N_clus = dim(pdf_array_1)[2]
  permn_list = combinat::permn(1:N_clus)
  
  min_dist = Inf
  for (permn in permn_list) {
    res = get_dist_betw_pdfarray(pdf_array_1[permn, permn, ], pdf_array_2, symmetric=TRUE, t_unit = t_unit)
    dist = res$dist
    n0_BetwSubj_mat = res$n0_mat
    
    if(dist < min_dist){
      the.permn = permn
      min_dist = dist
      the.n0_BetwSubj_mat = n0_BetwSubj_mat
    }
  }
  
  return(list(permn = the.permn, dist = min_dist, n0_mat = the.n0_BetwSubj_mat))
}
