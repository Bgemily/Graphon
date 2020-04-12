

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
  
  
  #k-means
  for (. in 1:MaxIter){
    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_mat=n0_mat, t_vec=t_vec, bw=bw, center_pdf_array = NA)
    clusters = res$clusters
    n0_mat = res$n0_mat
  }
  
  # get center_pdf_array
  n0_vec = get_n0_vec(n0_mat, clusters)
  center_pdf_array = get_center_pdf_array(edge_time_mat, clusters, n0_vec)
  
  
  return(list(clusters=clusters, n0_mat=n0_mat, center_pdf_array=center_pdf_array, clusters_specc=clusters_specc, n0_mat_specc=n0_mat_specc, clusters_specc_exaclus=clusters_specc_exaclus))
  
}

