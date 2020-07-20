

# main algorithm
do_cluster = function(edge_time_mat, N_clus, N_overclus=N_clus, MaxIter=10, N_trial=10, t_vec=seq(0, 50, 0.05), bw=1){
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(edge_time_mat)
  
  # initialize clusters via kmeanspp, initialize n0_mat 
  node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                      n0_vec = numeric(N_node), n0_mat = matrix(0,N_node,N_node), t_vec = t_vec, bw = bw)
  n0_vec_init = est_n0_vec(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), t_vec = t_vec, bw = bw)
  aligned_pdf_mat = t(sapply(1:N_node, function(i)shift(node_pdf_array[i,1,], n0_vec_init[i], pp=TRUE)))
  
  ## exact-clustering
  # res_exaclus = kmeanspp(data_mat = aligned_pdf_mat, N_clus = N_clus, MaxIter = MaxIter, N_trial = N_trial)
  # clusters_history = list(res_exaclus$clusters)

  ## over-clustering
  res_overclus = kmeanspp(data_mat = aligned_pdf_mat, N_clus = N_overclus, MaxIter = MaxIter, N_trial = N_trial)
  clusters = res_overclus$clusters
  res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  n0_mat = res$n0_mat
  n0_vec = res$n0_vec
  

  clusters_history = c(clusters_history, list(clusters))
  
  # merge clusters
  clusters_merged = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                   N_clus = N_clus, t_vec = t_vec, bw = bw)
  clusters_history = c(clusters_history, list(clusters_merged))
  
  
  
  # continue over-clustering
  for (. in 1:3){
    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec, n0_mat=n0_mat,
                         t_vec=t_vec, bw=bw, center_pdf_array = NULL)
    clusters = res$clusters
    n0_vec = res$n0_vec
    n0_mat = res$n0_mat

    clusters_history = c(clusters_history, list(clusters))
    
    # merge clusters
    clusters_merged = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                     N_clus = N_clus, t_vec = t_vec, bw = bw)
    clusters_history = c(clusters_history, list(clusters_merged))
    
  }
  
  
  # merge clusters
  clusters = clusters_merged
  res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  n0_mat = res$n0_mat
  n0_vec = res$n0_vec
  

  # continue k-means
  clusters_old = NULL
  n_iter = 1
  
  while (!identical(clusters, clusters_old) & n_iter<=MaxIter){
    clusters_old = clusters

    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec, n0_mat=n0_mat,
                         t_vec=t_vec, bw=bw, center_pdf_array = NULL)
    clusters = res$clusters
    n0_vec = res$n0_vec
    n0_mat = res$n0_mat

    n_iter = n_iter+1
    clusters_history = c(clusters_history, list(clusters))
  }


  
  # # or stick with spectral clustering
  # clusters_history = list(clusters)
  # for (. in 1:MaxIter){
  #   res = cluster_specc(edge_time_mat, clusters, n0_vec, N_clus, t_vec, bw)
  #   clusters = res$clusters
  #   n0_vec = res$n0_vec
  #   clusters_history = c(clusters_history, list(clusters))
  # }
  
  
  # get center_pdf_array
  center_pdf_array = get_center_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                          n0_vec = n0_vec, n0_mat = n0_mat, t_vec = t_vec, bw = bw)
  
  
  return(list(clusters=clusters, clusters_history=clusters_history, n0_vec=n0_vec, n0_mat=n0_mat, center_pdf_array=center_pdf_array))
  
}

