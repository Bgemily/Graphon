plot_pdf_array = function(pdf_array)
{
  k1 = dim(pdf_array)[1]
  k2 = dim(pdf_array)[2]
  par(mfrow = c(k1,k2))
  for (i in 1:k1) {
    for (j in 1:k2) {
      pdf = pdf_array[i,j,]
      plot(t, tail(pdf, length(t)), type='l', ylim = c(0,0.35))
    }
  }
  
  par(mfrow = c(1,1)); par(mfrow = c(1,1))
  # return()
}
