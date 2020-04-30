

# main algorithm
do_cluster = function(edge_time_mat, N_clus, N_overclus=N_clus, MaxIter=10, N_trial=10, t_vec=seq(0, 50, 0.05), bw=1){
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(edge_time_mat)
  
  # # init by k-medoids
  # N_node = nrow(edge_time_mat)
  # node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
  #                                     n0_vec = numeric(N_node), t_vec = t_vec, bw = bw)
  # n0_vec_init = est_n0_vec(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), t_vec = t_vec, bw = bw)
  # aligned_pdf_mat = t(sapply(1:N_node, function(i)shift(node_pdf_array[i,1,], n0_vec_init[i], pp=TRUE)))
  # dist_mat = rdist::pdist(aligned_pdf_mat)
  # clusters = mem2clus(cluster::pam(x = dist_mat, k = N_overclus, diss=TRUE)$cluster)
  # n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  # res_overclus=NULL
  # 
  # # cluster by k-medoids
  # clusters_history = list(clusters)
  # clusters_old = NULL
  # n_iter=0
  # while (!identical(clusters, clusters_old) & n_iter<=MaxIter){
  #   node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
  #                                        n0_vec = n0_vec, t_vec = t_vec, bw = bw)
  #   degree_mat = get_node_degree_mat(edge_time_mat = edge_time_mat, clusters = clusters)
  #   dist_mat = pairwise_dist_mat(node_pdf_array, degree_mat = degree_mat)$dist_mat
  #   
  #   clusters = mem2clus(cluster::pam(x = dist_mat, k = N_clus, diss=TRUE)$cluster)
  #   n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  #   
  #   n_iter = n_iter+1
  #   clusters_history = c(clusters_history, list(clusters))
  # }
  # 
  


  # # over-clustering via kmeanspp and no merge
  # clusters = list(c(1:nrow(edge_time_mat)))
  # n0_vec = numeric(nrow(edge_time_mat))
  # 
  # res_overclus = cluster_kmeans_overclus(edge_time_mat, clusters, n0_vec, N_overclus, N_clus = N_overclus,
  #                                       MaxIter = MaxIter, N_trial = N_trial, t_vec=t_vec, bw=bw)
  # 
  # clusters = res_overclus$clusters
  # n0_vec = res_overclus$n0_vec
  
  
  
  # naive initialization via kmeanspp
  node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                      n0_vec = numeric(N_node), t_vec = t_vec, bw = bw)
  n0_vec_init = est_n0_vec(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), t_vec = t_vec, bw = bw)
  aligned_pdf_mat = t(sapply(1:N_node, function(i)shift(node_pdf_array[i,1,], n0_vec_init[i], pp=TRUE)))
  
  # exact-clustering
  res_exaclus = kmeanspp(data_mat = aligned_pdf_mat, N_clus = N_clus, MaxIter = MaxIter, N_trial = N_trial)
  clusters_history = list(res_exaclus$clusters)
  
  # over-clustering
  res_overclus = kmeanspp(data_mat = aligned_pdf_mat, N_clus = N_overclus, MaxIter = MaxIter, N_trial = N_trial)
  clusters = res_overclus$clusters
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)

  clusters_history = c(clusters_history, list(clusters))
  
  # merge clusters
  clusters_merged = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, 
                                   N_clus = N_clus, t_vec = t_vec, bw = bw)
  clusters_history = c(clusters_history, list(clusters_merged))
  
  
  
  # continue over-clustering
  for (. in 1:3){
    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec,
                         t_vec=t_vec, bw=bw, center_pdf_array = NULL)
    clusters = res$clusters
    n0_vec = res$n0_vec

    clusters_history = c(clusters_history, list(clusters))
    
    # merge clusters
    clusters_merged = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, 
                                     N_clus = N_clus, t_vec = t_vec, bw = bw)
    clusters_history = c(clusters_history, list(clusters_merged))
    
  }
  
  
  # merge clusters
  clusters = clusters_merged
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  
  

  # continue k-means
  clusters_old = NULL
  n_iter = 1
  
  while (!identical(clusters, clusters_old) & n_iter<=MaxIter){
    clusters_old = clusters

    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec,
                         t_vec=t_vec, bw=bw, center_pdf_array = NULL)
    clusters = res$clusters
    n0_vec = res$n0_vec

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
                                          n0_vec = n0_vec, t_vec = t_vec, bw = 0.4)
  
  
  return(list(clusters=clusters, clusters_history=clusters_history, n0_vec=n0_vec, center_pdf_array=center_pdf_array))
  
}

