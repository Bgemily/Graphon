
plot_pdf_array = function(pdf_array, pdf_true_array = NULL, pdf_array_2 = NULL, t_vec = seq(0, 50, 0.05))
{
  k1 = dim(pdf_array)[1]
  k2 = dim(pdf_array)[2]
  par(mfrow = c(k1,k2))
  
  if(!is.null(pdf_true_array) & k1==k2){
    if (dim(pdf_array)[1]!=dim(pdf_true_array)[1] | dim(pdf_array)[2]!=dim(pdf_true_array)[2])
      stop("size of pdf_array and pdf_true_array should match.")
    
    permn = find_permn(pdf_true_array, pdf_array)$permn # find permutation for pdf_array
    pdf_array = pdf_array[permn, permn, ]
  }
  
  
  for (q in 1:k1) {
    for (l in 1:k2) {
      pdf_est = pdf_array[q,l,]
      
      if (length(pdf_est)!=length(t_vec))
        stop("Length of pdf and that of t_vec should match.")
      
      plot(t_vec, pdf_est, type='l', ylim = c(0,0.4))
      
      if(!is.null(pdf_array_2))
        lines(t_vec, pdf_array_2[q,l,], ylim = c(0, 0.4))
      
      if(!is.null(pdf_true_array)){
        lines(t_vec, pdf_true_array[q,l,], col='red', ylim=c(0,0.4))
      }
    }
  }
  
  par(mfrow = c(1,1)); par(mfrow = c(1,1))
  # return()
}
