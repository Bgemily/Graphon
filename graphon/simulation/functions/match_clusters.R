
# match clusters so that *clusters are of the same order* across subjects
match_clusters = function(edge_time_mat_list=NULL, clusters_list, n0_vec_list=NULL, n0_mat_list=NULL, pdf_array_list=NULL, pdf_true_array=NULL, t_vec=seq(0, 50, 0.05), bw=1){
  N_subj = length(clusters_list)
  t_unit = t_vec[2] - t_vec[1]
  
  if (N_subj==1) return(clusters_list)
  
  if (is.null(pdf_array_list))
    pdf_array_list = mapply(get_center_pdf_array, edge_time_mat = edge_time_mat_list, clusters = clusters_list, 
                            n0_vec = n0_vec_list, n0_mat = n0_mat_list, t_vec = list(t_vec), bw = list(bw), SIMPLIFY = FALSE)  
  if (is.null(pdf_true_array)){
    # warning: have not consider the case that pdf_array_list[[1]] has *less N_clus* than others, 
    # or pdf_array_list[[1]] is not a good estimate of pdf_true_array
    pdf_true_array = pdf_array_list[[1]] # use the first pdf_array as standard 
  }
  
  
  # Make all pdf_array have the same size as pdf_true_array
  subj_id = which(lapply(pdf_array_list, function(x)dim(x)[1]) < dim(pdf_true_array)[1])
  if (length(subj_id)>0){
    for (i in subj_id) {
      pdf_array = pdf_array_list[[i]]
      tmp = array(0, dim(pdf_true_array))
      tmp[1:dim(pdf_array)[1], 1:dim(pdf_array)[2], ] = pdf_array
      pdf_array_list[[i]] = tmp
    }
  }
  
  
  # find permutations, and align pdf_array across subjects
  permn_list = list() 
  for (s in 1:N_subj) {
    res = find_permn(pdf_array_1=pdf_array_list[[s]], pdf_array_2=pdf_true_array, t_unit = t_unit, t_vec=t_vec)
    permn = res$permn
    clusters_list[[s]] = clusters_list[[s]][permn]
    permn_list[[s]] = permn
    pdf_array_list[[s]] = pdf_array_list[[s]][permn, permn, ]
    
    n0_mat = res$n0_mat
    index_mat = matrix(1:length(n0_mat), dim(n0_mat))
    tmp.pdf.array = apply(index_mat, c(1,2), function(i){index = which(index_mat==i, arr.ind = TRUE); shift(f_origin = pdf_array_list[[s]][index[1], index[2], ], 
                                                                   n0 = n0_mat[index[1], index[2]], pp = TRUE)})
    pdf_array_list[[s]] = aperm(tmp.pdf.array, perm = c(2:ncol(n0_mat),1))
  }
  
  return(list(clusters_list=clusters_list, permn_list=permn_list, pdf_array_list=pdf_array_list))
}
