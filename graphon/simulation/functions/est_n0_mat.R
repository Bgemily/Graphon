# update n0_mat
# structure of n0_mat: tau_ij = vi*L_ql + vj*L_lq, where v is n0_vec, L is leading relation matrix

est_n0_mat = function(edge_time_mat, clusters, t_vec=seq(0,50,0.05), bw=1){
  
  N_clus = length(clusters)
  N_node = nrow(edge_time_mat)
  time_unit = t_vec[2]-t_vec[1]
  
  # estimate n0_vec
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  
  # estimate L_mat and n0_mat
  L_mat = matrix(data = NA, nrow = N_clus, ncol = N_clus)
  n0_mat = matrix(data = NA, nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      edge_time_submat = edge_time_mat[clusters[[q]], clusters[[l]], drop=F]
      
      
      tau_q_vec = time_unit * n0_vec[clusters[[q]]]
      tau_l_vec = time_unit * n0_vec[clusters[[l]]]
      
      
      edge_time_submat_1 = sweep(edge_time_submat, 1, tau_q_vec) # shift according to cluster q
      var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
      
      edge_time_submat_2 = sweep(edge_time_submat, 2, tau_l_vec) # shift according to cluster l
      var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
      
      
      if(var_1<var_2||is.na(var_1)) {
        # edge_time_submat = edge_time_submat_1
        n0_mat[clusters[[q]], clusters[[l]]] = matrix(n0_vec[clusters[[q]]], nrow = length(clusters[[q]]), ncol = length(clusters[[l]]))
        n0_mat[clusters[[l]], clusters[[q]]] = t(n0_mat[clusters[[q]], clusters[[l]]])
        L_mat[q,l] = 1; L_mat[l,q] = 0; if(q==l) L_mat[q,l] = 1
      }
      else {
        # edge_time_submat = edge_time_submat_2
        n0_mat[clusters[[q]], clusters[[l]]] = matrix(n0_vec[clusters[[l]]], nrow = length(clusters[[q]]), ncol = length(clusters[[l]]), byrow = TRUE)
        n0_mat[clusters[[l]], clusters[[q]]] = t(n0_mat[clusters[[q]], clusters[[l]]])
        L_mat[q,l] = 0; L_mat[l,q] = 1; if(q==l) L_mat[q,l] = 1
      }
      
      if(q==l){
        for (i in clusters[[q]]) {
          for (j in clusters[[l]]) {
            n0_mat[i,j] = min(n0_vec[i], n0_vec[j])
          }
        }
      }
      
    }
  }
  
  
  return(list(n0_mat=n0_mat, L_mat=L_mat, n0_vec=n0_vec))
  
}

