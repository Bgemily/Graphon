
combn_subj_recluster = function(subj_list, group_size = 5, MaxIter=2)
{
  for (group in split(1:length(subj_list), 1:(length(subj_list)/group_size))) {
    
    for (i in group) subj_list[[i]]$clusters_old = subj_list[[i]]$clusters
    
    for(. in 1:MaxIter){
      res = combn_subj(subj_list[group])
      pdf_array = res$pdf_array_combn
      for (i in group) {
        tmp = subj_list[[i]]
        permn = res$permn_list[[which(group==i)]]
        edge_time_mat = tmp$network$edge_time_mat; clusters = tmp$clusters; n0_vec=tmp$n0_vec
        
        subj_list[[i]]$clusters = re_cluster_connpatt(edge_time_mat, clusters, n0_vec, permn, pdf_array)$clusters
      }
    }
  }
  return(subj_list)
}

