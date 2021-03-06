
# Obtain a pdf from edge_time_vec. No time shifts.
get_pdf_vec = function(edge_time_vec, t_vec=seq(0, 50, 0.05), bw=1, intensity=TRUE){
  edge_time_vec = as.vector(edge_time_vec)
  
  jumps = edge_time_vec[which(edge_time_vec<=max(t_vec))]
  
  if (length(jumps) == 0)
    f_density = function(x)return(0)
  else{
    if (length(jumps) == 1) 
      jumps = c(jumps, jumps + 1e-10)
    d = density(jumps, bw = bw)
    f_density = approxfun(d, yleft = 0, yright = 0)
  }
  pdf_vec = sapply(t_vec, f_density)
  
  
  if(intensity) {
    conn_prob = sum(edge_time_vec<=max(t_vec))/length(edge_time_vec)
    pdf_vec = pdf_vec * conn_prob
  }
  
  
  return(pdf_vec)
}

