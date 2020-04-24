# update n0_vec
est_n0_vec = function(edge_time_mat, clusters, t_vec=seq(0,50,0.05), bw=1){
  
  t_unit = t_vec[2] - t_vec[1]
  
  # do not shift any event at begining
  n0_vec = numeric(nrow(edge_time_mat))
  
  node_pdf_array = get_node_pdf_array(edge_time_mat, clusters, n0_vec, t_vec, bw)
  for (l in 1:length(clusters)) {
    if (length(clusters[[l]]) == 1) {
      n0_vec[clusters[[l]]] = 0
      next
    }
    
    # align each node with the first node (i_0) in this cluster
    i_0 = clusters[[l]][1]
    for (i in clusters[[l]]) {
      # align i-th row towards i_0-th row
      possible_n0 = get_dist_betw_pdfarray(node_pdf_array[i, , , drop=F], node_pdf_array[i_0, , , drop=F], symmetric=FALSE, t_unit = t_unit)$n0_mat
      if(length(clusters)>1) possible_n0 = possible_n0[-l] # only consider time shifts in connection with other clusters. This is crucial.
      n0_vec[i] = possible_n0[which.max(abs(possible_n0))]
    }
    # make the minimum time shift to be zero
    n0_vec[clusters[[l]]] = n0_vec[clusters[[l]]] - min(n0_vec[clusters[[l]]])
  }
  
  return(n0_vec)
}