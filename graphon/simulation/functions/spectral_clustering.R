
pairwise_dist_mat = function(f_list, pp=FALSE, n0_mat=matrix(0,length(f_list),length(f_list))){
  dist_mat = matrix(0, length(f_list), length(f_list))
  for (i in 1:(length(f_list)-1)) {
    for (j in ((i+1):length(f_list))) {
      res = align_curves_gd(f_list[[i]], f_list[[j]], n0=n0_mat[i,j], step_size = 0.02, pp=pp)
      dist_mat[i,j] = res$dist_min
      dist_mat[j,i] = dist_mat[i,j]
      
      n0_mat[i,j] = res$n0
      n0_mat[j,i] = -n0_mat[i,j]
    }
  }
  return(list(dist_mat=dist_mat, n0_mat=n0_mat))
}

pairwise_corr_mat = function(f_list, pp=FALSE, n0_mat=matrix(0,length(f_list),length(f_list))){
  res = pairwise_dist_mat(f_list, pp, n0_mat)
  dist_mat = sqrt(res$dist_mat)
  corr_mat = exp(-dist_mat^2/median(dist_mat)^2)
  return(list(corr_mat = corr_mat, n0_mat = res$n0_mat, dist_mat = dist_mat))
}

spectral_clustering = function (W, k) 
{
  n = ncol(W)
  S = rowSums(W)
  D = diag(S)
  # L = D - W
  D_sq_inv = diag(S^(-1/2))
  # L = (diag(1, n) - D_sq_inv%*%W%*%D_sq_inv)
  L = (diag(1, n) - diag(S^(-1))%*%W)
  U = (eigen(L)$vectors)[, ((n - k + 1):n)]
  C = cluster::pam(x = U, k = k)
  return(C$clustering)
}