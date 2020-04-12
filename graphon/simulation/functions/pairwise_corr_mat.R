

pairwise_corr_mat = function(pdf_array, degree_mat){
  res = pairwise_dist_mat(pdf_array, degree_mat)
  dist_mat = sqrt(res$dist_mat)
  corr_mat = exp(-dist_mat^2/median(dist_mat)^2)
  return(list(corr_mat = corr_mat, n0_mat = res$n0_mat, dist_mat = dist_mat))
}

