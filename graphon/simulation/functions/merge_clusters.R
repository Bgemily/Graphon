
# obtain connecting pattern (pdf of shifted events) for each pair of clusters
get_pdf_array_ = function(clusters, edge_time_mat, n0_vec, t=seq(0, 50, 0.05), h=1)
{  
  time_unit = t[2]-t[1]
  pdf_array = array(dim=c(length(clusters),length(clusters),length(t)))
  for (i in 1:length(clusters)) {
    for (j in 1:length(clusters)) {
      edge_time_submat = edge_time_mat[clusters[[i]], clusters[[j]]]
      
      if (sum(edge_time_submat<Inf)==0) {
        pdf = rep(0, 2*length(t))
      }
      else 
      {
        tau_i_vec = time_unit * n0_vec[clusters[[i]]]
        tau_j_vec = time_unit * n0_vec[clusters[[j]]]
        
        tau_ij_mat_1 = matrix(tau_i_vec, length(tau_i_vec), length(tau_j_vec)) # shift according to cluster i
        edge_time_submat_1 = matrix(edge_time_submat-tau_ij_mat_1, nrow=1)
        var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
        tau_ij_mat_2 = t(matrix(tau_j_vec, length(tau_j_vec), length(tau_i_vec))) # shift according to cluster j
        edge_time_submat_2 = matrix(edge_time_submat-tau_ij_mat_2, nrow=1)
        var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
        
        if((var_1<var_2) || is.na(var_1)){
          tau_ij_mat = tau_ij_mat_1
          edge_time_submat = edge_time_submat_1
        }
        else{
          tau_ij_mat = tau_ij_mat_2
          edge_time_submat = edge_time_submat_2
        }
        
        # tau_ij_mat = matrix(tau_i_vec, length(tau_i_vec), length(tau_j_vec))
        # edge_time_submat = matrix(edge_time_submat-tau_ij_mat, nrow=1)
        
        pdf = get_pp_f_list(edge_time_submat, t, h=h)[[1]]
      }
      
      pdf_array[i,j,] = tail(pdf, length(t))
    }
  }
  
  return(pdf_array)
}


# obtain distance between each pair of clusters
get_dist_mat_ = function(pdf_array)
{
  if (dim(pdf_array)[1]!=dim(pdf_array)[2]) {
    stop("The dim of pdf_array should be k*k*n!")
  }
  dist_mat_  = matrix(Inf,dim(pdf_array)[1], dim(pdf_array)[1]) #Diagonal is Inf because we need min non-zero distance when merge.
  for (i in 1:(dim(pdf_array)[1]-1)) {
    for (j in (i+1):dim(pdf_array)[2]) {
      
      dist_mat_[i,j] = sum(sapply(c(1:dim(pdf_array)[2]), function(l)align_pdf_gd(pdf_array[i,l,], pdf_array[j,l,])$dist_min))
      dist_mat_[j,i] = dist_mat_[i,j]
      # 
      # f_row_i = matrix(t(pdf_array[i,,]), nrow=1)
      # f_row_j = matrix(t(pdf_array[j,,]), nrow=1)
      # dist_mat_[i,j] = norm(f_row_j-f_row_i, type="2")
      # 
      # f_col_i = matrix(t(pdf_array[,i,]), nrow=1)
      # f_col_j = matrix(t(pdf_array[,j,]), nrow=1)
      # dist_mat_[j,i] = norm(f_col_j-f_col_i, type="2")
    }
  }
  dist_mat_ = 1/2*(t(dist_mat_)+dist_mat_)
  return(dist_mat_)
}


merge_clusters = function(clusters, edge_time_mat, n0_mat, k, t=seq(0,50,0.05), h=1)
{
  n0_vec = rep(0, length(f_list))
  for (i in 1:length(clusters)) {
    n0_vec[clusters[[i]]] = n0_mat[clusters[[i]], clusters[[i]][1]] - min(n0_mat[clusters[[i]], clusters[[i]][1]])
  }
  pdf_array = get_pdf_array_(clusters, edge_time_mat, n0_vec, t, h)
  while (length(clusters)>k) {
    dist_mat_ = get_dist_mat_(pdf_array)
    merge_clus_ind = which(dist_mat_==min(dist_mat_), arr.ind = TRUE)[1,]
    
    clusters[[merge_clus_ind[1]]] = c(clusters[[merge_clus_ind[1]]], clusters[[merge_clus_ind[2]]])
    clusters[[merge_clus_ind[2]]] = NULL
    
    n0_vec = rep(0, length(f_list))
    for (i in 1:length(clusters)) {
      n0_vec[clusters[[i]]] = n0_mat[clusters[[i]], clusters[[i]][1]] - min(n0_mat[clusters[[i]], clusters[[i]][1]])
    }
    pdf_array = get_pdf_array_(clusters, edge_time_mat, n0_vec, t, h)
  }
  return(clusters)
}