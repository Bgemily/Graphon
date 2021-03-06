### Initialization of clusters and n0_mat
### ERROR: this function should take ALL subjects as input.
get_init_v2 = function(edge_time_mat_list, N_clus, t_vec=seq(0, 50, 0.05), order_true_list=NULL){
  
  time_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  # Initialize time shifts -------------------------------------------------------

  n0_vec_list = list()
  aligned_cdf_longmat = c()
  for (m in 1:N_subj) {
    edge_time_mat = edge_time_mat_list[[m]]
    N_node = nrow(edge_time_mat)
    node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                           n0_mat = matrix(0,N_node,N_node), t_vec = t_vec)
    earliest_edge_time = apply(edge_time_mat, 1, function(row) min(row[which(row>1)]))
    
    n0_vec = (earliest_edge_time)/time_unit
    n0_vec = round(n0_vec)
    
    n0_vec = n0_vec - min(n0_vec)
    # n0_mat = n0_vec2mat(n0_vec = n0_vec)
    
    ### Debug
    if (!is.null(order_true_list)) {
      n0_vec[order_true_list[[m]]] = sort(n0_vec)
    }
    
    
    aligned_cdf_mat = t(sapply(1:N_node, function(i) shift_v2(f_origin=node_cdf_array[i,1,], n0=-n0_vec[i])))
    
    n0_vec_list[[m]] = n0_vec
    aligned_cdf_longmat = rbind(aligned_cdf_longmat, aligned_cdf_mat)
  }
  

  # Initialize clusters -----------------------------------------------------

  
  res = cluster::pam(x=aligned_cdf_longmat, k=N_clus, diss=FALSE, cluster.only=TRUE)
  
  membership_list = list()
  clusters_list = list()
  for (m in 1:N_subj) {
    row_id = 1:N_node_vec[m] + I(m>=2)*sum(N_node_vec[1:(m-1)])
    membership = res[row_id]
    clusters = mem2clus(membership)
    
    membership_list[[m]] = membership
    clusters_list[[m]] = clusters
  }
  
  
  
  return(list(membership_list=membership_list, 
              clusters_list=clusters_list, 
              n0_vec_list=n0_vec_list))
}

### Test
# edge_time_mat = matrix(1:4,2,2)
# edge_time_mat = kronecker(edge_time_mat, matrix(10,5,5))
# edge_time_mat_list = list(edge_time_mat, edge_time_mat)
# res = get_init_v2(edge_time_mat_list = edge_time_mat_list, N_clus = 2)
# res$clusters_list
# res$n0_vec_list
# res$membership_list
