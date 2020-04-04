plot_pdf_array = function(pdf_array)
{
  k = dim(pdf_array)[1]
  par(mfrow = c(k,k))
  for (i in 1:k) {
    for (j in 1:k) {
      pdf = pdf_array[i,j,]
      plot(t, tail(pdf, length(t)), type='l', ylim = c(0,0.35))
    }
  }
  
  par(mfrow = c(1,1)); par(mfrow = c(1,1))
  # return()
}
