
# match clusters so that *clusters are of the same order* across subjects
match_clusters = function(edge_time_mat_list, clusters_list, n0_vec_list, t_vec=seq(0, 50, 0.05), bw=1){
  N_subj = length(clusters_list)
  
  if (N_subj==1) return(clusters_list)
  
  pdf_array_list = mapply(get_center_pdf_array, edge_time_mat_list, clusters_list, n0_vec_list, list(t_vec), list(bw), SIMPLIFY = FALSE)  
  
  for (s in 2:N_subj) {
    permn = find_permn(pdf_array_1=pdf_array_list[[1]], pdf_array_2=pdf_array_list[[s]])$permn
    clusters_list[[s]] = clusters_list[[s]][permn]
  }
  
  return(clusters_list)
}
