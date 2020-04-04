obtain_pdf = function(cdf_vec, t, bw=1){
  if (length(cdf_vec)!=length(t)) stop("Length of cdf and t do not match!")
  
  jumps = t[which(cdf_vec[-1]-cdf_vec[-length(cdf_vec)] > 0)]
  if (length(jumps)==1) jumps = c(jumps, jumps+1e-10)
  d = density(jumps, bw=bw)
  smooth.pdf = approxfun(d, yleft=0, yright=0)
  return(list(density = d, jumps=jumps, smooth.pdf=smooth.pdf))
}