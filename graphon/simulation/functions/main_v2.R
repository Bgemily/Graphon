

main_v2 = function(SEED=NULL, N_subj=1, N_node_vec = rep(90,N_subj),
                   N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
                   total_time=200, 
                   conn_patt_var=1, conn_patt_sep = 1.5, const=40, conn_prob_mean = 1, conn_prob_rad = 0, 
                   time_shift_struc=max, time_shift_mean_vec = rep(10,N_clus), time_shift_rad = 10,
                   ...)
{
  ### Generate networks
  network_list = generate_network2_v2(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                             N_clus = N_clus, clus_size_mat = clus_size_mat,
                             total_time = total_time, 
                             conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                             conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                             time_shift_struc = time_shift_struc,
                             time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)

  edge_time_mat_list = network_list$edge_time_mat_list
  cdf_true_array = network_list$cdf_true_array
  pdf_true_array = network_list$pdf_true_array
  membership_true_list = network_list$membership_true_list
  v_true_list = network_list$time_shift_list

  
  ### Apply algorithm
  res = do_cluster_v2(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus, 
                      total_time = total_time, max_iter=1, ...)
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$center_cdf_array -> center_cdf_array_est
  
  
  ### Compute accuracy of clusters, time shifts and conn patt
  ARI_vec = numeric(length = N_subj)
  for (m in 1:N_subj) {
    ARI_tmp = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est[[m]]), 
                          memb_true_vec = membership_true_list[[m]])
    ARI_vec[m] = ARI_tmp
  }
  ARI_mean = mean(ARI_vec)
  
  
  normed_dist_mat = matrix(nrow=N_clus, ncol=N_clus)
  unnormed_dist_mat = matrix(nrow=N_clus, ncol=N_clus)
  conn_prob_err_mat = matrix(nrow=N_clus, ncol=N_clus)
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      conn_prob_est = tail(center_cdf_array_est[q,k,], 1)
      conn_prob_true = tail(cdf_true_array[q,k,], 1)
      
      cdf_est_normed = center_cdf_array_est[q,k,] / conn_prob_est
      cdf_true_normed = cdf_true_array[q,k,] / conn_prob_true
      
      normed_dist_mat[q,k] = sqrt(sum( (cdf_est_normed - cdf_true_normed)^2 ))
      unnormed_dist_mat[q,k] = sqrt(sum( (center_cdf_array_est[q,k,] - cdf_true_array[q,k,])^2 ))
      conn_prob_err_mat[q,k] = abs(conn_prob_est - conn_prob_true)
    }
  }
  
  
  spearman_corr_vec = numeric(length = N_subj)
  pearson_corr_vec = numeric(length = N_subj)
  for (m in 1:N_subj) {
    spearman_corr_vec[m] = cor(v_true_list[[m]], v_vec_list_est[[m]], method = "spearman")
    pearson_corr_vec[m] = cor(v_true_list[[m]], v_vec_list_est[[m]], method = "pearson")
  }
  spearman_corr_vec_mean = mean(spearman_corr_vec)
  pearson_corr_vec_mean = mean(pearson_corr_vec)

    
  
  # return(list(network=network, clus_result=res))
  network_param = list(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                       N_clus = N_clus, clus_size_mat = clus_size_mat,
                       total_time = total_time, 
                       conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                       conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                       time_shift_struc = time_shift_struc,
                       time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  
  return(list(network_param=network_param, 
              ARI_vec=ARI_vec, ARI_mean=ARI_mean,
              normed_dist_mat=normed_dist_mat, unnormed_dist_mat=unnormed_dist_mat, conn_prob_err_mat=conn_prob_err_mat,
              spearman_corr_vec=spearman_corr_vec, spearman_corr_vec_mean=spearman_corr_vec_mean,
              pearson_corr_vec=pearson_corr_vec, pearson_corr_vec_mean=pearson_corr_vec_mean))
}




