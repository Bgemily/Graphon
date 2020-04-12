

# k-means clustering
cluster_kmeans = function(edge_time_mat, clusters, n0_mat, center_pdf_array=NA, t_vec=seq(0, 50, 0.05), bw=1){
  
  N_node = nrow(edge_time_mat); N_clus = length(clusters)
  
  n0_vec = get_n0_vec(n0_mat, clusters)
  node_pdf_array = get_node_pdf_array(edge_time_mat, clusters, n0_vec, t_vec, bw)
  
  if (is.na(center_pdf_array)) {
    center_pdf_array = get_center_pdf_array(edge_time_mat, clusters, n0_vec, t_vec, bw)
  }
  
  degree_mat = get_node_degree_mat(edge_time_mat, clusters)
  membership = numeric(N_node)
  for (i in 1:N_node) {
    if(sum(degree_mat[i,])>0) weights = degree_mat[i,]/sum(degree_mat[i,])
    else weights = NA
    
    dist_vec = numeric(N_clus)
    for (l in 1:N_clus) {
      # distance between ith node and l-th cluster
      dist_vec[l] = get_dist_betw_pdfarray(node_pdf_array[i, , , drop=F], center_pdf_array[l, , ,drop=F], symmetric=FALSE, weights=weights)$dist
    }
    membership[i] = which.min(dist_vec)
  }
  
  clusters = mem2clus(membership)
  
  return(list(clusters = clusters, n0_mat=n0_mat))
}

