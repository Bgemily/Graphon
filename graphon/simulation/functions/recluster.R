# recluster nodes based on current pdf_array(:the estimated connecting pattern matrix) and clusters and permutation
re_cluster_connpatt = function(edge_time_mat, clusters, n0_vec, permn, pdf_array, time_unit=0.05, t=seq(0, 50, 0.05), h=1){
  clusters = clusters[permn]
  clusters_upd = vector('list', length(clusters))
  
  
  tmp = function(i){
    pdf_array_node_i = pdf_array[1, , ,drop=F]
    
    for (l in 1:dim(pdf_array)[2]) {
      edge_time_submat = edge_time_mat[i,clusters[[l]], drop=FALSE]
      
      edge_time_submat_1 = edge_time_submat
      var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
      
      tau_vec = time_unit * n0_vec[clusters[[l]]]
      edge_time_submat_2 = edge_time_submat-tau_vec
      var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
      
      if(var_1<var_2 || is.na(var_1)){
        pdf = get_pp_f_list(edge_time_submat_1, t=t, h=h)[[1]]
        pdf_array_node_i[1, l, ] = tail(pdf, length(t))
      }
      else{
        pdf = get_pp_f_list(edge_time_submat_2, t=t, h=h)[[1]]
        pdf_array_node_i[1, l, ] = tail(pdf, length(t))
      }
      
    }
    
    dist_vec = rep(Inf, dim(pdf_array)[1])
    for (l in 1:dim(pdf_array)[1]) {
      # dist: distance between ith node and l-th cluster
      dist = get_dist_betw_PatternMat_(pdf_array_node_i, pdf_array[l, , ,drop=F], symmetric=FALSE)$dist
      dist_vec[l] = dist
    }

    return(which.min(dist_vec))
  }
  
  
  membership = sapply(1:nrow(edge_time_mat), tmp)
  
  for (l in 1:length(clusters)) {
    clusters_upd[[l]] = which(membership==l)
  }

      
  return(list(clusters = clusters_upd))
}