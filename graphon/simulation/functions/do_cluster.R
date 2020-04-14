

# main algorithm
do_cluster = function(edge_time_mat, N_clus, N_overclus=N_clus, MaxIter=3, t_vec=seq(0, 50, 0.05), bw=1){
  
  # initial clustering (spectral clustering)
  clusters = list(c(1:nrow(edge_time_mat)))
  n0_mat = matrix(0, nrow(edge_time_mat),  ncol(edge_time_mat))
  
  # res = cluster_specc(edge_time_mat, clusters, n0_mat, N_clus, t_vec, bw)
  res = cluster_specc_overclus(edge_time_mat, clusters, n0_mat, N_overclus, N_clus, t_vec, bw)
  
  clusters = res$clusters
  n0_mat = res$n0_mat
  
  # save the clustering result of spectral clustering
  clusters_specc = res$clusters
  n0_mat_specc = res$n0_mat
  clusters_specc_exaclus = res$clusters_exaclus
  
  
  # k-means
  for (. in 1:MaxIter){
    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_mat=n0_mat, t_vec=t_vec, bw=bw, center_pdf_array = NA)
    clusters = res$clusters
    n0_mat = res$n0_mat
  }
  
  
  # or stick with spectral clustering
  # for (. in 1:MaxIter){
  #   res = cluster_specc(edge_time_mat, clusters, n0_mat, N_clus, t_vec, bw)
  #   clusters = res$clusters
  #   n0_mat = res$n0_mat
  # }
  
  
  # update n0_vec
  n0_vec = numeric(nrow(n0_mat))
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
      possible_n0 = get_dist_betw_pdfarray(node_pdf_array[i, , , drop=F], node_pdf_array[i_0, , , drop=F], symmetric=FALSE)$n0_mat
      possible_n0 = possible_n0[-l] # only consider time shifts in connection with other clusters. This is crucial.
      n0_vec[i] = possible_n0[which.max(abs(possible_n0))]
    }
    # make the minimum time shift to be zero
    n0_vec[clusters[[l]]] = n0_vec[clusters[[l]]] - min(n0_vec[clusters[[l]]])
  }
  
  # get center_pdf_array
  center_pdf_array = get_center_pdf_array(edge_time_mat, clusters, n0_vec)
  
  
  return(list(clusters=clusters, n0_mat=n0_mat, n0_vec=n0_vec, center_pdf_array=center_pdf_array, clusters_specc=clusters_specc, n0_mat_specc=n0_mat_specc, clusters_specc_exaclus=clusters_specc_exaclus))
  
}

