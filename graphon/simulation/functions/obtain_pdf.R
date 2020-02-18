obtain_pdf = function(cdf_vec, t){
  if (length(cdf_vec)!=length(t)) stop("Length of cdf and t do not match!")
  
  jumps = t[which(cdf_vec[-1]-cdf_vec[-length(cdf_vec)] > 0)]
  d = density(jumps, bw='SJ')
  return(list(density = d, jumps=jumps))
}