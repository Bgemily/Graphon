

# k-means clustering
cluster_kmeans = function(edge_time_mat, clusters, n0_vec, n0_mat=NULL, center_pdf_array=NULL, t_vec=seq(0, 50, 0.05), bw=1, intensity=TRUE){
  
  n0_mat = NULL ###### change this to use new method
  
  N_node = nrow(edge_time_mat); 
  N_clus = length(clusters)
  if(!is.null(center_pdf_array) && length(clusters)!=dim(center_pdf_array)[2]) 
    stop("dim(center_pdf_array)[2] and length(clusters) should match.")
  else if(!is.null(center_pdf_array) && length(clusters)!=dim(center_pdf_array)[1])
    N_clus = dim(center_pdf_array)[1]
  
  t_unit = t_vec[2]-t_vec[1]
  
  # estimate node_pdf_array
  node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                      t_vec = t_vec, bw = bw)
  
  # estimate center_pdf_array
  if (is.null(center_pdf_array)) {
    center_pdf_array = get_center_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                            t_vec = t_vec, bw = bw)
  }
  
  # update membership
  degree_mat = get_node_degree_mat(edge_time_mat = edge_time_mat, clusters = clusters, intensity=intensity) 
  membership = numeric(N_node)
  dist_to_centr_vec = numeric(N_node)
  for (i in 1:N_node) {
    if(sum(degree_mat[i,])>0) {
      # NEED JUSTIFICATION
      # weights = log(degree_mat[i,]+1)
      weights = sqrt(degree_mat[i,]) ############
      weights = weights/sum(weights)
    }
    else weights = NULL
    
    dist_vec = numeric(N_clus)
    for (l in 1:N_clus) {
      # distance between ith node and l-th cluster
      dist_vec[l] = get_dist_betw_pdfarray(pdf_array_1 = node_pdf_array[i, , , drop=F], 
                                           pdf_array_2 = center_pdf_array[l, , ,drop=F], 
                                           symmetric=FALSE, weights=weights, t_unit=t_unit, t_vec = t_vec)$dist
    }
    membership[i] = which.min(dist_vec)
    dist_to_centr_vec[i] = min(dist_vec)
  }
  
  clusters = mem2clus(membership)
  
  # update n0_vec
  # n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  n0_mat = res$n0_mat
  n0_vec = res$n0_vec
  
  return(list(clusters = clusters, n0_vec=n0_vec, n0_mat=n0_mat, dist_to_centr_vec=dist_to_centr_vec))
}

