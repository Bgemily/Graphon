
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

distance2 = function(theta_prime, gamma_prime, n0) # get rid of high frequency,same result as distance() 
{
  N = length(theta_prime)
  k = seq(1,N-1)
  k_freq = k
  k_freq[(N%/%2):(N-1)] = k_freq[(N%/%2):(N-1)]-N
  return(sum(abs(theta_prime[k+1]*exp(1i*2*pi*k_freq*n0/N)-gamma_prime[k+1])^2)+abs(theta_prime[1]+n0-gamma_prime[1])^2)
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


# Align curve (using grid search) -----------------------------------------

align_curves_search_ = function(f_origin, f_shift, n_min, n_max, by=1)
{
  if (n_max < n_min) stop("n_max shoule be no less than n_min.")
  if (n_min < 0) stop("n_min should be negative.")
  theta = fft(f_origin)
  gamma = fft(f_shift)
  N = length(theta)
  k = seq(1, N-1)
  theta_prime = c(theta[1], theta[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  gamma_prime = c(gamma[1], gamma[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  
  n0 = n_min
  dist_min = Inf
  for (n in seq(n_min, n_max, by)) {
    d = distance(theta_prime, gamma_prime, n)
    if (d<dist_min){
      n0 = n
      dist_min = d
    }
  }
  return(list(n0=n0, dist_min=dist_min))
}


# Allows to shift towards left or right, both n_min and n_max can be negative.
align_curves_search = function(f1, f2, n_min, n_max, by=1)
{
  if (n_max >= n_min && n_min >= 0)
    return(align_curves_search_(f1, f2, n_min, n_max, by))
  else if (n_min <= n_max && n_max <= 0)
  {
    r = align_curves_search_(f2, f1, -n_max, -n_min, by)
    return(list(n0 = -r$n0, dist_min = r$dist_min))
  }
  else if (n_min <= 0 && 0 <= n_max)
  {
    r1 = align_curves_search_(f1, f2, 0, n_max, by)
    r2 = align_curves_search_(f2, f1, 0, -n_min, by)
    if (r1$dist_min <= r2$dist_min)
      return(r1)
    else
      return(list(n0=-r2$n0, dist_min=r2$dist_min))
  }
}

# Shift the grid search region if the result is near n_min or n_max
align_curves_search2 = function(f1, f2, n_min, n_max, by=100)
{
  r = align_curves_search(f1, f2, n_min, n_max, by)
  if (r$n0-by >= n_min && r$n0+by <= n_max)
    return(r)
  else if (r$n0-by < n_min)
    return(align_curves_search2(f1, f2, 2*n_min-n_max, n_min+2*by, by))
  else
    return(align_curves_search2(f1, f2, n_max-2*by, n_max*2-n_min, by))
}

# align_curves_search(f_origin, f_shift, 0, 2000)
# align_curves_search(f_shift, f_origin, -300,300)


# Align curve (using gradient) -------------------------------------------------------------

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

# Allow n0 to be negative
align_curves_gd = function(f1, f2, n0, step_size, MaxIter=1000, stopping_redu=0.01, pp=FALSE)
{
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


# align_curves_gd(f_origin, f_origin, -2000, .1, stopping_redu = 0.01)


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


# Find first time f>0+ and last time f<1- -----------------------------------

find_mint_ = function(f, min_thres = 0.05)
{
  value = f[1]; i = 1
  while (value<min_thres) {
    i = i+1
    value = f[i]
  }
  return(i)
}

find_maxt_ = function(f, max_thres=0.95)
{
  i = length(f); value = f[i]
  while (value > max_thres) {
    i = i-1
    value = f[i]
  }
  return(i)
}

# find_maxt_(f_origin)
# find_mint_(f_origin)





# # Test --------------------------------------------------------------------
# 
# x = seq(0,35,0.01)
# f_origin = c(rep(0,length(x)),x/max(x))
# N = length(f_origin)
# plot(f_origin, type = 'l')
# 
# theta = fft(f_origin)
# 
# tmp = c(rep(0,length(theta))); tmp[2:3]=theta[2:3]
# plot(Re(fft(tmp, inverse = T)),type = 'l', ylim = c(-2000,3000))
# par(new=T)
# tmp2 = c(rep(0,length(theta))); tmp2[4]=theta[4]
# plot(Re(fft(tmp2, inverse = T)),type = 'l',ylim = c(-2000,3000), ylab='')
# par(new=T)
# plot(Re(fft(tmp+tmp2, inverse = T)),type = 'l', ylim = c(-2000,3000),col=2, ylab='')
# 
# 
# plot(abs(theta))
# 
# f_shift = shift(f_origin, 200)
# # f_shift = c(rep(0,length(x)-200),rep(1/length(x), length(x)), rep(0,length(x)+200))
# plot(f_shift, type='l')
# 
# theta_shift = fft(f_shift)
# 
# 
# ##### Test coef formula
# # {
# #   k = seq(0,N-1)
# #   n0 = 2/0.01
# #   theta_hope = theta*exp(1i*2*pi*k*(2/0.01)/N)+exp(-1i*2*pi*k)*(exp(1i*2*pi*n0*k/N)-1)/(1-exp(-1i*2*pi*k/N))
# #   theta_hope[1] = theta[1]+n0
# #   
# #   plot(abs(theta_hope-theta_shift))
# #   diff = theta_hope-theta_shift
# #   diff[1]
# # }
# 
# 
# ##### plot of dist between curves
# 
# {
#   theta = fft(f_origin)
#   gamma = fft(f_shift)
#   N = length(theta)
#   k = seq(1, N-1)
#   theta_prime = c(theta[1], theta[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
#   gamma_prime = c(gamma[1], gamma[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
#   
#   seq = seq(0,400,.1)
#   dists_vec = sapply(seq, function(x)distance2(theta_prime, gamma_prime, x))
#   plot(c(seq),c(dists_vec), type = 'l', xlab = 'integer n (time shift)', ylab = 'distance between two curves')
#   abline(h=distance(theta_prime, gamma_prime, 200), col=2)
# 
# }
# 
