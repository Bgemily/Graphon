
# pdf_array:  n*k*N, degree_mat: n*k
pairwise_dist_mat = function(pdf_array, degree_mat){
  N_node = dim(pdf_array)[1]
  dist_mat = matrix(0, N_node, N_node)
  n0_mat = matrix(0, N_node, N_node)
  for (i in 1:(N_node-1)) {
    for (j in ((i+1):N_node)) {
      weights = degree_mat[i,]*degree_mat[j,]; weights = weights/sum(weights)
      
      res = get_dist_betw_pdfarray(pdf_array[i, , , drop=F], pdf_array[j, , , drop=F], symmetric=FALSE, weights=weights)
      # res = align_curves_gd(f_list[[i]], f_list[[j]], n0=n0_mat[i,j], step_size = 0.02, pp=pp)
      
      dist_mat[i,j] = res$dist
      dist_mat[j,i] = dist_mat[i,j]
      
      n0 = max(res$n0_mat) ### need justification
      n0_mat[i,j] = n0
      n0_mat[j,i] = -n0_mat[i,j]
    }
  }
  return(list(dist_mat=dist_mat, n0_mat=n0_mat))
}

