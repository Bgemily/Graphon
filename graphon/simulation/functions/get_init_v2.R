### Initialization of clusters and n0_mat
get_init_v2 = function(edge_time_mat, N_clus, t_vec=seq(0, 50, 0.05)){
  
  time_unit = t_vec[2] - t_vec[1]
  N_node = nrow(edge_time_mat)
  
# Initialize time shifts -------------------------------------------------------

  node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                      n0_mat = matrix(0,N_node,N_node), t_vec = t_vec)
  earliest_edge_time = apply(edge_time_mat, 1, function(row) min(row[which(row>1)]))
  
  n0_vec = (earliest_edge_time)/time_unit
  n0_vec = round(n0_vec)
  
  n0_vec = n0_vec - min(n0_vec)
  # n0_mat = n0_vec2mat(n0_vec = n0_vec)

# Initialize clusters -----------------------------------------------------

  aligned_cdf_mat = t(sapply(1:N_node, function(i) shift_v2(f_origin=node_cdf_array[i,1,], n0=-n0_vec[i])))
  
  res = cluster::pam(x=aligned_cdf_mat, k=N_clus, diss=FALSE, cluster.only=TRUE)
  clusters = mem2clus(res)
  
  
  return(list(membership=res, clusters=clusters, n0_vec=n0_vec))
}

### Test
# edge_time_mat = matrix(1:4,2,2)
# edge_time_mat = kronecker(edge_time_mat, matrix(10,5,5))
# res = get_init_v2(edge_time_mat = edge_time_mat, N_clus = 2)
# res$clusters
# res$n0_vec
# res$membership
