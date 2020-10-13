eval_loss = function(edge_time_mat, n0_mat, clusters, center_cdf_array=NULL, 
                     standardize=FALSE, t_vec=seq(0,50,0.05)){

  N_node = nrow(edge_time_mat)
  N_clus = length(clusters)
  t_unit = t_vec[2]-t_vec[1]
  node_cdf_array = get_node_cdf_array(edge_time_mat = edge_time_mat, clusters = mem2clus(1:N_node),
                                      n0_mat = n0_mat, t_vec = t_vec, standardize = standardize)
  if(is.null(center_cdf_array)){
    center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                            n0_mat = n0_mat, t_vec = t_vec, standardize = standardize)
  }
  
  loss = 0
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      tmp = sweep(node_cdf_array[clusters[[q]],clusters[[l]],], MARGIN = 3, STATS = center_cdf_array[q,l,])
      loss = loss + sum(tmp^2)*t_unit
    }
  }
  loss = loss/N_node^2/max(t_vec) # loss: mean squared error, invariant to N_Node and T.
  
  return(list(loss=loss))
}


# res = results1$conn_prob_0.2[[1]]
# eval_loss(edge_time_mat = res$network$edge_time_mat,n0_mat = res$clus_result$n0_mat,
#           clusters=res$clus_result$clusters,t_vec = res$network$t_vec,standardize = F)
# eval_loss(edge_time_mat = res$network$edge_time_mat,n0_mat = res$network$tau_mat/res$network$t_vec[2],
#           clusters=mem2clus(res$network$membership_true),t_vec = res$network$t_vec,standardize = F)

