
# Shift curve towards left -----------------------------------------------

shift = function(f_origin, n0, pad = 1, pp=FALSE)
{
  N = length(f_origin)
  if (!pp)
    return(c(f_origin[(n0+1):N], rep(pad, n0)))
  else
    return(c(f_origin[(n0+1):N], rep(0, n0)))
}



# Compute distance with given shift ---------------------------------------

distance = function(theta_prime, gamma_prime, n0, pp=FALSE) # squared L2 distance
{
  N = length(theta_prime)
  k = seq(1,N-1)
  if (!pp)
    return(sum(abs(theta_prime[k+1]*exp(1i*2*pi*k*n0/N)-gamma_prime[k+1])^2)+abs(theta_prime[1]+n0-gamma_prime[1])^2)
  else
    return(sum(abs(theta_prime[k+1]*exp(1i*2*pi*k*n0/N)-gamma_prime[k+1])^2))
}



# gradient -----------------------------------------------------------------

gradient = function(theta_prime, gamma_prime, n0, pp=FALSE)
{
  N = length(theta_prime)
  if(N != length(gamma_prime)) stop("Length of theta and gamma do not match.")
  k = seq(1, N-1)
  k_freq = k
  k_freq[(N%/%2):(N-1)] = k_freq[(N%/%2):(N-1)]-N
  
  if (!pp)
    return ((4*pi/N) * sum(k_freq*Im((theta_prime[k+1])*Conj(gamma_prime[k+1])*exp(1i*2*pi*k_freq*n0/N))) + 2*n0 + 2*Re(theta_prime[1]-gamma_prime[1]) )
  else
    return ((4*pi/N) * sum(k_freq*Im((theta_prime[k+1])*Conj(gamma_prime[k+1])*exp(1i*2*pi*k_freq*n0/N))) )
}

# grads = sapply(seq(0,400,0.1),function(x)gradient(theta_prime, gamma_prime, x))
# plot(seq(0,400,0.1), grads, type = 'l')
# abline(h=0,col=2)


# Align curve (using gradient) -------------------------------------------------------------

# Non-negative n0. Move f_origin towards f_shift by n0. If pp=FALSE, pad 1 at the end of f_origin when shifting; otherwise, pad 0.
align_curves_gd_ = function(f_origin, f_shift, n0, step_size, MaxIter=1000, stopping_redu=0.01, pp=FALSE)
{
  theta = fft(f_origin)
  gamma = fft(f_shift)
  N = length(theta)
  k = seq(1, N-1)
  {
    if (!pp){
    theta_prime = c(theta[1], theta[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
    gamma_prime = c(gamma[1], gamma[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
    }
    else{
      theta_prime = theta
      gamma_prime = gamma
    }
  }
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  while (dist_redu>stopping_redu && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd = gradient(theta_prime, gamma_prime, n0, pp=pp); 
    n0 = n0 - step_size*gd
    n0 = round(n0) 
    
    if(n0<0) { 
      n0 = 0; 
      dist_curr = distance(theta_prime, gamma_prime, n0, pp=pp)
      break 
    }
    if(n0>length(f_origin)) { 
      n0 = length(f_origin); 
      dist_curr = distance(theta_prime, gamma_prime, n0, pp=pp)
      break 
    }
    
    dist_upd = distance(theta_prime, gamma_prime, n0, pp=pp)
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    dist_curr = dist_upd
  }
  
  # print(iter_count)
  
  return(list(n0=n0, dist_min = dist_curr))
}

# Allow n0 to be negative. Move f1 towards f2 by n0. Pad 0 at the head of f1 and f2. 
align_curves_gd = function(f1, f2, n0, step_size, MaxIter=1000, stopping_redu=0.01, pp=FALSE)
{
  f1 = c(rep(0, length(f1)), f1)
  f2 = c(rep(0, length(f2)), f2)

  
  if (n0 <= 0) {
    r = align_curves_gd_(f2, f1, -n0, step_size, MaxIter, stopping_redu, pp=pp)
    if(r$n0>0)  {return(list(n0 = -r$n0, dist_min = r$dist_min))}
    else {return(align_curves_gd_(f1, f2, 0, step_size, MaxIter, stopping_redu, pp=pp))}
  }
  else{
    r = align_curves_gd_(f1, f2, n0, step_size, MaxIter, stopping_redu, pp=pp)
    if(r$n0>0) {return(r)}
    else {
      r = (align_curves_gd_(f2, f1, 0, step_size, MaxIter, stopping_redu, pp=pp))
      return(list(n0 = -r$n0, dist_min = r$dist_min))
    }
  }
}


# Initialize n0 by aligning cdf1 and cdf2.
align_pdf_gd = function(pdf1, pdf2, n0=0, step_size=0.02, MaxIter=1000, stopping_redu=0.01)
{
  if(sum(pdf1)==0 || sum(pdf2)==0) n0_init=n0
  else{
    cdf1 = cumsum(pdf1)/sum(pdf1)
    cdf2 = cumsum(pdf2)/sum(pdf2)
    n0_init = align_curves_gd(cdf1, cdf2, n0, step_size, MaxIter, stopping_redu, pp=FALSE)$n0
  }
  res = align_curves_gd(pdf1, pdf2, n0_init, step_size, MaxIter, stopping_redu, pp=TRUE)
  return(list(n0 = res$n0, dist_min = res$dist_min))
}


