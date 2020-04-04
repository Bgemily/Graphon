get_dist_betw_PatternMat_ = function(pdf_array_1, pdf_array_2, symmetric=TRUE)
{ 
  if(symmetric){
    ### assuming the connecting pattern matrix is symmetric and square.
    dist = 0
    n0_BetwSubj_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    for (i in 1:(dim(pdf_array_1)[1])) {
      for (j in (i):dim(pdf_array_1)[2]) {
        res = align_pdf_gd(pdf_array_1[i,j,], pdf_array_2[i,j,])
        
        dist = dist + res$dist_min
        n0_BetwSubj_mat[i,j] = res$n0
        n0_BetwSubj_mat[j,i] = n0_BetwSubj_mat[i,j]
      }
    }
    return(list(dist = dist, n0_BetwSubj_mat = n0_BetwSubj_mat))
  }
  
  else{
    dist = 0
    n0_BetwSubj_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    for (i in 1:(dim(pdf_array_1)[1])) {
      for (j in (1):dim(pdf_array_1)[2]) {
        res = align_pdf_gd(pdf_array_1[i,j,], pdf_array_2[i,j,])
        
        dist = dist + res$dist_min
        n0_BetwSubj_mat[i,j] = res$n0
      }
    }
    return(list(dist = dist, n0_BetwSubj_mat = n0_BetwSubj_mat))
  }
    
}

match_clus = function(pdf_array_1, pdf_array_2) # find the optimal permutation
{
  if (dim(pdf_array_1)[1]!=dim(pdf_array_2)[1]) {
    stop("Shapes of pdf_array_1 and pdf_array_2 should be the same!")
  }
  k = dim(pdf_array_1)[1]
  permn_list = combinat::permn(1:k)
  
  min.dist = Inf
  for (permn in permn_list) {
    res = get_dist_betw_PatternMat_(pdf_array_1, pdf_array_2[permn, permn, ])
    dist = res$dist
    n0_BetwSubj_mat = res$n0_BetwSubj_mat
    
    if(dist < min.dist){
      the.permn = permn
      min.dist = dist
      the.n0_BetwSubj_mat = n0_BetwSubj_mat
    }
  }
  
  return(list(the.permn = the.permn, min.dist = min.dist, the.n0_BetwSubj_mat = the.n0_BetwSubj_mat))
}


# simply combine nodes in the same cluster of two (matched) subjects
merge_clus_ = function(clusters1, clusters2){
  if (length(clusters1)!=length(clusters2)) {
    stop("Two clustering should have the same number of clusters!")
  }
  
  # Add nodes of subject2 to clusters1 by making their index range from (1+#nodes of subj1) to (#nodes of subj2+#nodes of subj1)
  clusters_combn = lapply(1:length(clusters1), function(i) c(clusters1[[i]], clusters2[[i]]+length(unlist(clusters1))))
  return(clusters_combn)
}




# combine clusters and edge_time_mat and n0_vec. clusters1 and clusters2 may have different number of nodes
combn_clus = function(clusters1, edge_time_mat1, n0_vec1, clusters2, edge_time_mat2, n0_vec2){
  pdf_array_1 = get_pdf_array_(clusters1, edge_time_mat1, n0_vec1)
  pdf_array_2 = get_pdf_array_(clusters2, edge_time_mat2, n0_vec2)
  
  
  permn = match_clus(pdf_array_1, pdf_array_2)$the.permn
  clusters_combn = merge_clus_(clusters1, clusters2[permn])
  
  edge_time_combn = matrix(Inf, nrow=nrow(edge_time_mat1)+nrow(edge_time_mat2), ncol=ncol(edge_time_mat1)+ncol(edge_time_mat2))
  edge_time_combn[1:nrow(edge_time_mat1), 1:ncol(edge_time_mat1)] = edge_time_mat1
  edge_time_combn[1:nrow(edge_time_mat2)+nrow(edge_time_mat1), 1:ncol(edge_time_mat2)+ncol(edge_time_mat1)] = edge_time_mat2
  
  # does not take into account the time lag difference between subjects
  # i.e. assume the distribution of time lags are consistent among subjects
  n0_vec_combn = c(n0_vec1,n0_vec2)
  
  
  return(list(clusters_combn=clusters_combn, edge_time_combn=edge_time_combn, n0_vec_combn=n0_vec_combn, permn = permn))
}


# combine subjects and update pdf_array
combn_subj = function(subj_list){
  if(length(subj_list)==0){
    stop("Empty subject list!")
  }
  
  s1 = subj_list[[1]]
  clusters1 = s1$clusters
  edge_time_mat1 = s1$network$edge_time_mat
  n0_vec1 = s1$n0_vec
  
  permn_list = vector("list", length=length(subj_list))
  permn_list[[1]] = (c(1,2,3))
  
  if(length(subj_list)==1){
    pdf_array_combn = get_pdf_array_(clusters1, edge_time_mat1, n0_vec1)
    return(list(pdf_array_combn = pdf_array_combn, permn_list=permn_list))
  }
  
  
  for (i in 2:length(subj_list)) {
    s2 = subj_list[[i]]
    clusters2 = s2$clusters
    edge_time_mat2 = s2$network$edge_time_mat; 
    n0_vec2 = s2$n0_vec
    
    res = combn_clus(clusters1, edge_time_mat1, n0_vec1, clusters2, edge_time_mat2, n0_vec2)
    clusters1 = res$clusters_combn
    edge_time_mat1 = res$edge_time_combn
    n0_vec1 = res$n0_vec_combn
    permn_list[[i]] = res$permn
  }
  

  pdf_array_combn = get_pdf_array_(clusters1, edge_time_mat1, n0_vec1)
  
  return(list(pdf_array_combn = pdf_array_combn, permn_list=permn_list))
}




