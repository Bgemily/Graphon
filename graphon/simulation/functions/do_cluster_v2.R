

# main algorithm
do_cluster_v2 = function(edge_time_mat_list, N_clus, 
                         total_time = 200, t_vec=seq(0,total_time,length.out=1000),
                         MaxIter=10, conv_thres=1e-3, ...)
{
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  clusters_history = list()
  cluster_time = align_time = 0
  
# Initialize clusters and time shifts -------------------------------------
  
  clusters_list = list()
  n0_vec_list = list()
  n0_mat_list = list()
  for (m in 1:N_subj) {
    edge_time_mat = edge_time_mat_list[[m]]
    
    ### Get initialization
    res = get_init_v2(edge_time_mat=edge_time_mat, N_clus=N_clus, t_vec=t_vec)
    clusters = res$clusters
    n0_vec = res$n0_vec
    n0_mat = n0_vec2mat(n0_vec = n0_vec)
    
    clusters_list[[m]] = clusters
    n0_vec_list[[m]] = n0_vec
    n0_mat_list[[m]] = n0_mat
    
  }
  

  clusters_history = c(clusters_history, list(clusters_list))
  
  
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, t_vec = t_vec)
  
  

# Update clusters and connecting patterns alternatively ---------------------

  n_iter = 1
  stopping = FALSE
  
  clusters_list_update = clusters_list_current = clusters_list
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  center_cdf_array_update = center_cdf_array_current = center_cdf_array
  
  while (!stopping & n_iter<=MaxIter){
    
    ### Update clusters, time shifts and connecting patterns
    res = cluster_kmeans_v2(edge_time_mat_list=edge_time_mat_list, 
                            clusters_list=clusters_list_current, 
                            n0_vec_list=n0_vec_list_current, n0_mat_list=n0_mat_list_current, 
                            center_cdf_array = center_cdf_array_current,
                            t_vec=t_vec, ...)
    clusters_list_update = res$clusters_list
    n0_vec_list_update = res$n0_vec_list
    n0_mat_list_update = res$n0_mat_list
    v_vec_list_update = res$v_vec_list
    center_cdf_array_update = res$center_cdf_array
    
    ### Record computing time for clustering and aligning
    cluster_time = cluster_time + res$cluster_time
    align_time = align_time + res$align_time
    
    
    clusters_history = c(clusters_history, list(clusters_list_update))
    
  
    ### Evaluate stopping criterion
    delta_n0_vec = sum((unlist(n0_vec_list_update)-unlist(n0_vec_list_current))^2) / 
                    ( sum(unlist(n0_vec_list_current)^2) + .Machine$double.eps )

    delta_clusters_vec = numeric(length = N_subj)
    for (m in 1:N_subj) {
      clusters_update = clusters_list_update[[m]]
      clusters_current = clusters_list_current[[m]]
      delta_tmp = 1 - get_one_ARI(memb_est_vec = clus2mem(clusters_update), 
                                       memb_true_vec = clus2mem(clusters_current))
      delta_clusters_vec[m] = delta_tmp
    }
    weights = N_node_vec/sum(N_node_vec)
    delta_clusters = sum(delta_clusters_vec * weights)

    delta_center_cdf = tryCatch(sum((center_cdf_array_update-center_cdf_array_current)^2) / 
                                  (sum(center_cdf_array_current^2) + .Machine$double.eps),
                                error=function(x)1)
    stopping = mean(c(delta_center_cdf,delta_clusters,delta_n0_vec)) < conv_thres
    
    
    ### *update -> *current
    n_iter = n_iter+1
    clusters_list_update -> clusters_list_current
    n0_vec_list_update -> n0_vec_list_current 
    n0_mat_list_update -> n0_mat_list_current 
    v_vec_list_update -> v_vec_list_current
    center_cdf_array_update -> center_cdf_array_current 
    
  }

  
  ### Get final result
  clusters_list_current -> clusters_list
  n0_vec_list_current -> n0_vec_list
  n0_mat_list_current -> n0_mat_list
  v_vec_list_current -> v_vec_list
  center_cdf_array_current -> center_cdf_array
  
  center_pdf_array = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, t_vec = t_vec)
  
  
  return(list(clusters_list=clusters_list, clusters_history=clusters_history, 
              v_vec_list=v_vec_list,
              center_pdf_array=center_pdf_array, center_cdf_array=center_cdf_array,
              cluster_time=cluster_time, align_time=align_time))
  
}


# Test --------------------------------------------------------------------

# res = generate_network2_v2(SEED=10, N_subj=2, N_node_vec=rep(90,2),
#                            N_clus=3,
#                            total_time=200,
#                            conn_patt_var=1, conn_patt_sep = 1.5, const=40, conn_prob_mean = 1, conn_prob_rad = 0,
#                            time_shift_struc=max, time_shift_rad = 10)
# 
# edge_time_mat_list = res$edge_time_mat_list
# cdf_true_array = res$cdf_true_array
# pdf_true_array = res$pdf_true_array
# clusters_list = res$clus_true_list
# time_shift_list = res$time_shift_list
# 
# tmp = do_cluster_v2(edge_time_mat_list = edge_time_mat_list, N_clus = 3, max_iter=1)
# tmp$clusters_list
# tmp$clusters_history
# plot(time_shift_list[[1]], tmp$v_vec_list[[1]])

# tmp2 = do_cluster(edge_time_mat = edge_time_mat_list[[1]], N_clus = 3)
# tmp2$clusters
# tmp2$clusters_history[[1]]
# plot(time_shift_list[[1]], tmp2$n0_vec)
