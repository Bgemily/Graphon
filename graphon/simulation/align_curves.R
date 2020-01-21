
# Shift curve towards right -----------------------------------------------

shift = function(f_origin, n0)
{
  N = length(f_origin)
  return(c(f_origin[(n0+1):N], rep(1,n0)))
}



# Compute distance with given shift ---------------------------------------

distance = function(theta_prime, gamma_prime, n0) # squared L2 distance
{
  N = length(theta_prime)
  k = seq(1,N-1)
  return(sum(abs(theta_prime[k+1]*exp(1i*2*pi*k*n0/N)-gamma_prime[k+1])^2)+abs(theta_prime[1]+n0-gamma_prime[1])^2)
}


# gradient -----------------------------------------------------------------

gradient = function(theta_prime, gamma_prime, n0)
{
  N = length(theta_prime)
  if(N != length(gamma_prime)) stop("Length of theta and gamma do not match.")
  k = seq(1, N-1)
  return (-(4*pi/N) * sum(k*Re(1i*(theta_prime[k+1])*Conj(gamma_prime[k+1])*exp(1i*2*pi*k*n0/N))) + 2*n0 + 2*Re(theta_prime[1]-gamma_prime[1]) )
}

grads = sapply(seq(199,201,0.1),function(x)gradient(theta_prime, gamma_prime, x))
plot(seq(199,201,0.1), grads, type = 'l')
abline(h=0,col=2)


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


# Align curve (using gradient) -------------------------------------------------------------

align_curves_gd = function(f_origin, f_shift, n0, step_size)
{
  theta = fft(f_origin)
  gamma = fft(f_shift)
  N = length(theta)
  k = seq(1, N-1)
  theta_prime = c(theta[1], theta[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  gamma_prime = c(gamma[1], gamma[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  MaxIter = 30
  for (tmp in 1:MaxIter) {
    gd = gradient(theta_prime, gamma_prime, n0)
    n0 = n0 - step_size*sign(gd)
    n0 = round(n0)
    if(n0<0) n0 = 0
  }
  return(n0)
}

# align_curves_gd(f_origin, f_shift, 201, 1)



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



# Compute the mean curve ------------------------------

mean_curve = function(f_list, n0_vec)
{
  shifted_f_mat = matrix(0, nrow=length(f_list), ncol=length(f_list[[1]]))
  for (i in (1:length(f_list))) {
    shifted_f_mat[i,] = shift(f_list[[i]], n0_vec[i])
  }
  return(colMeans(shifted_f_mat))
}



# # Test --------------------------------------------------------------------
# 
# x = seq(0,35,0.01)
# f_origin = c(rep(0,length(x)),x/max(x))
# # f_origin = c(rep(0,length(x)),1/3*pnorm(x,5,2)+2/3*pnorm(x,10,3))
# N = length(f_origin)
# plot(f_origin)
# 
# theta = fft(f_origin)
# 
# f_shift = shift(f_origin, 200)
# plot(f_shift)
# 
# theta_shift = fft(f_shift)
# 
# 
# ##### Test coef formula
# 
# k = seq(0,N-1)
# n0 = 2/0.01
# theta_hope = theta*exp(1i*2*pi*k*(2/0.01)/N)+exp(-1i*2*pi*k)*(exp(1i*2*pi*n0*k/N)-1)/(1-exp(-1i*2*pi*k/N))
# 
# 
# plot(abs(theta_hope-theta_shift))
# diff = theta_hope-theta_shift
# diff[1]
# 
# ##### plot of dist between curves
# 
# theta = fft(f_origin)
# gamma = fft(f_shift)
# N = length(theta)
# k = seq(1, N-1)
# theta_prime = c(theta[1], theta[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
# gamma_prime = c(gamma[1], gamma[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
# 
# 
# seq = seq(0,400,1)
# # seq = c(seq,seq+.1,seq+.2,seq+.3,seq+.4,seq+.5,seq+.6,seq+.7,seq+.8,seq+.9)
# dists_vec = sapply(seq, function(x)distance(theta_prime, gamma_prime, x))
# dists_vec2 = sapply(seq2, function(x)distance(theta_prime, gamma_prime, x))
# plot(c(seq),c(dists_vec), type = 'l', xlab = 'integer n (time shift)', ylab = 'distance between two curves')
# abline(h=distance(theta_prime, gamma_prime, 200), col=2)
# 
# 
# ##### align curves
# 
# # align_curves_search(f_origin, f_shift, 0, 2000)
# # align_curves_search(f_shift, f_origin, -300,300)
# 
