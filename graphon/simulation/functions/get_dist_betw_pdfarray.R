
# compute distance between two pdf_array
get_dist_betw_pdfarray = function(pdf_array_1, pdf_array_2, symmetric=FALSE, weights=NULL, t_unit=0.05){ 
  if(symmetric){
    ### assuming the connecting pattern matrix is symmetric and square. Used for brute-force search of permutation.
    dist = 0
    n0_BetwSubj_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    for (i in 1:(dim(pdf_array_1)[1])) {
      for (j in (i):dim(pdf_array_1)[2]) {
        res = align_pdf_gd(pdf_array_1[i,j,], pdf_array_2[i,j,], t_unit=t_unit)
        
        dist = dist + res$dist_min
        n0_BetwSubj_mat[i,j] = res$n0
        n0_BetwSubj_mat[j,i] = n0_BetwSubj_mat[i,j]
      }
    }
    return(list(dist = dist, n0_mat = n0_BetwSubj_mat))
  }
  
  else{
    dist_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    n0_BetwSubj_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    for (i in 1:(dim(pdf_array_1)[1])) {
      for (j in (1):dim(pdf_array_1)[2]) {
        res = align_pdf_gd(pdf_array_1[i,j,], pdf_array_2[i,j,], t_unit=t_unit)
        
        dist_mat[i,j] = res$dist_min
        n0_BetwSubj_mat[i,j] = res$n0
      }
    }
    
    if(is.null(weights)) dist = mean(dist_mat)
    else dist = sum(dist_mat*weights)

    return(list(dist = dist, n0_mat = n0_BetwSubj_mat))
  }
  
}

