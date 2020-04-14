

# obtain connecting pattern (pdf of shifted events) for each pair of clusters
get_center_pdf_array = function(edge_time_mat, clusters, n0_vec, t_vec=seq(0, 50, 0.05), bw=1){  
  time_unit = t_vec[2]-t_vec[1]
  pdf_array = array(dim=c(length(clusters),length(clusters),length(t_vec)))
  for (i in 1:length(clusters)) {
    for (j in 1:length(clusters)) {
      edge_time_submat = edge_time_mat[clusters[[i]], clusters[[j]], drop=F]
      
      
      tau_i_vec = time_unit * n0_vec[clusters[[i]]]
      tau_j_vec = time_unit * n0_vec[clusters[[j]]]
      
      
      edge_time_submat_1 = sweep(edge_time_submat, 1, tau_i_vec) # shift according to cluster i
      var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
      
      edge_time_submat_2 = sweep(edge_time_submat, 2, tau_j_vec) # shift according to cluster j
      var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
      
      
      if(var_1<var_2||is.na(var_1)) edge_time_submat = edge_time_submat_1
      else edge_time_submat = edge_time_submat_2
      
      # NEED JUSTIFICATION
      # deal with the situation that shifted event times are negative (thus will be eliminated in estimated pdf)
      if (min(edge_time_submat) < min(t_vec)) {
        MAX_edge_time_submat = max(edge_time_submat[is.finite(edge_time_submat)])
        if (MAX_edge_time_submat-min(edge_time_submat) <= max(t_vec)-min(t_vec)) 
          edge_time_submat = edge_time_submat - min(edge_time_submat) + min(t_vec)
        else{ 
          N_early_events = sum( edge_time_submat < MAX_edge_time_submat-(max(t_vec)-min(t_vec)) )
          N_late_events = sum( (edge_time_submat > min(edge_time_submat)+max(t_vec)-min(t_vec)) & is.finite(edge_time_submat))
          if (N_early_events > N_late_events){ # keep lower tail
            edge_time_submat = edge_time_submat - min(edge_time_submat) + min(t_vec)
            edge_time_submat[edge_time_submat>max(t_vec)] = Inf
          }
          else{
            edge_time_submat = edge_time_submat + max(t_vec) - MAX_edge_time_submat
            edge_time_submat[edge_time_submat<0] = Inf
          }
        }
      }
      
      
      pdf_array[i,j,] = get_pdf_vec(edge_time_vec = edge_time_submat, t_vec, bw=bw)
    }
  }
  
  return(pdf_array)
}

