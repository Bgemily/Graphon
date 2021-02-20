# Turn n0_vec to n0_mat through maximum operations

n0_vec2mat = function(n0_vec){
  n0_mat = matrix(nrow=length(n0_vec), ncol=length(n0_vec))
  for (i in 1:nrow(n0_mat)) {
    for (j in 1:ncol(n0_mat)) {
      n0_mat[i,j] = max(n0_vec[i], n0_vec[j])
    }
  }
  return(n0_mat)
}