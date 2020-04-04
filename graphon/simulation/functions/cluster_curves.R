# source('~/Documents/Academic/SC/graphon/simulation/align_curves.R')
# source('./align_curves.R')


# Compute the mean curve ------------------------------

mean_curve = function(f_list, n0_vec, pp=FALSE)
{
  shifted_f_mat = matrix(0, nrow=length(f_list), ncol=length(f_list[[1]]))
  for (i in (1:length(f_list))) {
    shifted_f_mat[i,] = shift(f_list[[i]], n0_vec[i], pp=pp)
  }
  return(colMeans(shifted_f_mat))
}



# Kmeans++ initialization -------------------------------------------------

init_kmeans = function(f_list, k, seed, step_size=0.01, pp=FALSE)
{
  index = c()
  set.seed(seed); seed = seed+10; index = c(index,sample(length(f_list), 1))
  f_center_list = f_list[index]
  for (i in 2:k) {
    r = re_cluster_gd(f_list, f_center_list, rep(0,length(f_list)), step_size, pp=pp)
    dist_min_vec = r$dist_min_vec
    set.seed(seed); seed = seed+100.2; 
    N = 5
    candidates = tail(order(dist_min_vec), N)
    index = c(index,sample(candidates, 1, prob=dist_min_vec[candidates]/sum(dist_min_vec[candidates])))
    f_center_list = f_list[index]
  }
  return(index)
}

# init_kmeans(f_list, 3, 23)

# Re-cluster (using grid search) --------------------------------------------------------------

re_cluster = function(f_list, f_center_list, mint_f_vec, maxt_f_vec)
{
  z = vector(mode = "numeric", length = length(f_list))
  clusters = vector(mode = "list", length = length(f_center_list))
  for (i in (1:length(f_list))) {
    f = f_list[[i]]
    dist_min = Inf
    mint_f = mint_f_vec[i]; maxt_f = maxt_f_vec[i]
    for (l in (1:length(f_center_list))) {
      f_center = f_center_list[[l]]
      
      mint_f_center = find_mint_(f_center); maxt_f_center = find_maxt_(f_center)
      n_min = min(mint_f - mint_f_center, maxt_f - maxt_f_center)
      n_max = max(mint_f - mint_f_center, maxt_f - maxt_f_center)
      
      d = align_curves_search2(f, f_center, n_min, n_max, by=100)$dist_min
      if (d < dist_min){
        dist_min = d
        z[i] = l
      }
    }
    clusters[[z[i]]] = c(clusters[[z[i]]], i)
  }
  clusters = clusters[!sapply(clusters, is.null)]
  return(clusters) 
}

# mint_f_vec = sapply(f_list, find_mint_)
# maxt_f_vec = sapply(f_list, find_maxt_)
# re_cluster(f_list, f_center_list, mint_f_vec, maxt_f_vec)



# Re-cluster (using GD) ---------------------------------------------------

re_cluster_gd = function(f_list, f_center_list, n0_vec, step_size, pp=FALSE)
{
  z = vector(mode = "numeric", length = length(f_list))
  clusters = vector(mode = "list", length = length(f_center_list))
  dist_min_vec = vector(mode="numeric", length = length(f_list))
  for (i in (1:length(f_list))) {
    f = f_list[[i]]
    n0 = n0_vec[i]
    dist_min = Inf
    for (l in (1:length(f_center_list))) {
      f_center = f_center_list[[l]]
      
      r = align_curves_gd(f, f_center, n0, step_size, pp=pp)
      d = r$dist_min
      if (d < dist_min){
        dist_min = d
        z[i] = l
        n0_vec[i] = r$n0
      }
    }
    dist_min_vec[i] = dist_min
    clusters[[z[i]]] = c(clusters[[z[i]]], i)
  }
  clusters = clusters[!sapply(clusters, is.null)]
  return(list(clusters=clusters, n0_vec=n0_vec, dist_min_vec=dist_min_vec)) 
}

# re_cluster_gd(f_list, f_center_list, rep(0,length(f_list)), 0.05)



# Re-center ---------------------------------------------------------------

re_center = function(f_list, clusters, mint_f_vec, maxt_f_vec)
{
  f_center_list = vector(mode = "list", length = length(clusters))
  dist_vec = vector("numeric", length(clusters))
  for (k in (1:length(clusters))) {
    f_clusterk_list = f_list[clusters[[k]]]
    mint_clusterk_vec = mint_f_vec[clusters[[k]]]
    maxt_clusterk_vec = maxt_f_vec[clusters[[k]]]
    
    r = re_center_(f_clusterk_list, mint_clusterk_vec, maxt_clusterk_vec)
    f_center_list[[k]] = r$f_center
    dist_vec[k] = r$dist
  }
  return(list(f_center_list=f_center_list, mean_dist=mean(dist_vec)))
}

re_center_ = function(f_clusterk_list, mint_clusterk_vec, maxt_clusterk_vec, MaxIter = 100)
{
  f_center = f_clusterk_list[[which.min(mint_clusterk_vec)]]
  dist_curr = Inf
  dist_reduc = Inf
  iter_count = 0
  while (dist_reduc>0.01 && iter_count <= MaxIter) {
    iter_count = iter_count + 1
    n0_vec = vector("numeric", length(f_clusterk_list))
    dist_vec = vector("numeric", length(f_clusterk_list))
    mint_f_center = find_mint_(f_center); maxt_f_center = find_maxt_(f_center)
    
    for (i in (1:length(f_clusterk_list))) {
      mint_f = mint_clusterk_vec[i]; maxt_f = maxt_clusterk_vec[i]
      n_min = min(mint_f - mint_f_center, maxt_f - maxt_f_center)
      n_max = max(mint_f - mint_f_center, maxt_f - maxt_f_center)
      
      r = align_curves_search2(f_clusterk_list[[i]], f_center, n_min, n_max, by=100)
      n0_vec[i] = r$n0; dist_vec[i] = r$dist_min
    }
    dist_upd = mean(dist_vec)
    dist_reduc = (dist_curr-dist_upd) / (dist_upd)
    if (is.na(dist_reduc)) dist_reduc = 0
    dist_curr = dist_upd
    
    n0_vec = n0_vec - min(n0_vec)
    f_center = mean_curve(f_clusterk_list, n0_vec)
  }
  # print(iter_count)
  # plot(f_center, type = 'l')
  return(list(f_center=f_center, dist=dist_curr))
}


# clusters = list(c(1:5,16:20), c(6:15))
# re_center(f_list, clusters, mint_f_vec, maxt_f_vec)$mean_dist



# Re-center (using GD) ----------------------------------------------------

re_center_gd = function(f_list, clusters, n0_vec, step_size, MaxIter=100, pp=FALSE)
{
  f_center_list = vector(mode = "list", length = length(clusters))
  dist_vec = vector("numeric", length(clusters))
  for (k in (1:length(clusters))) {
    f_clusterk_list = f_list[clusters[[k]]]
    n0_clusterk_vec = n0_vec[clusters[[k]]]
    
    r = re_center_gd_(f_clusterk_list, n0_clusterk_vec, step_size, MaxIter, pp=pp)
    f_center_list[[k]] = r$f_center
    dist_vec[k] = r$dist
    n0_vec[clusters[[k]]] = r$n0_clusterk_vec
  }
  # mean_dist: mean with-in cluster distance
  return(list(f_center_list=f_center_list, n0_vec=n0_vec, mean_dist=mean(dist_vec)))
}



re_center_gd_ = function(f_clusterk_list, n0_clusterk_vec, step_size, MaxIter = 100, pp=FALSE)
{
  f_center = f_clusterk_list[[which.min(n0_clusterk_vec)]]
  dist_curr = Inf
  dist_reduc = Inf
  iter_count = 0
  while (dist_reduc>0.01 && iter_count <= MaxIter) {
    iter_count = iter_count + 1
    dist_vec = vector("numeric", length(f_clusterk_list))

    for (i in (1:length(f_clusterk_list))) { # align all f to the current f_center
      r = align_curves_gd(f_clusterk_list[[i]], f_center, n0=n0_clusterk_vec[i], step_size, pp=pp)
      n0_clusterk_vec[i] = r$n0; dist_vec[i] = r$dist_min
    }
    dist_upd = mean(dist_vec)
    dist_reduc = (dist_curr-dist_upd) / (dist_upd)
    if (is.na(dist_reduc)) dist_reduc = 0
    dist_curr = dist_upd
    
    n0_clusterk_vec = n0_clusterk_vec-min(n0_clusterk_vec)
    f_center = mean_curve(f_clusterk_list, n0_clusterk_vec, pp=pp)
  }
  # print(iter_count)
  # plot(f_center, type = 'l')
  return(list(f_center=f_center, dist=dist_curr, n0_clusterk_vec=n0_clusterk_vec))
}


# clusters = list(c(1:5,16:20), c(6:15))
# re_center_gd(f_list, clusters, rep(0,length(f_list)), 0.05)$mean_dist




# Cluster curves ----------------------------------------------------------

cluster_curves = function(f_list, k, MaxIter=100, seed=45, stopping_redu=0.01)
{
  mint_f_vec = sapply(f_list, find_mint_)
  maxt_f_vec = sapply(f_list, find_maxt_)
  set.seed(seed);
  index = sample(length(f_list), k)
  # index = c(1,2,3,4)
  f_center_list = f_list[index]
  
  dist_curr = Inf
  dist_redu = Inf
  iter_count = 0
  while(dist_redu > stopping_redu && iter_count <= MaxIter)
  {
    iter_count = iter_count+1
    
    clusters = re_cluster(f_list, f_center_list, mint_f_vec, maxt_f_vec)
    r = re_center(f_list, clusters, mint_f_vec, maxt_f_vec)
    f_center_list = r$f_center_list
    dist_upd = r$mean_dist
    
    dist_redu = (dist_curr-dist_upd) / (dist_upd)
    if (is.na(dist_redu)) dist_redu=0
    dist_curr = dist_upd
  }
  clusters = re_cluster(f_list, f_center_list, mint_f_vec, maxt_f_vec)
  # print(iter_count)
  # print(dist_curr)
  plot(f_center_list[[1]], type='l')
  plot(f_center_list[[2]], type='l')
  f_center_var = 0
  for (i in (1:(length(f_center_list)-1))) {
    for (j in ((i+1):length(f_center_list))) {
      f_center_var = f_center_var+align_curves_search2(f_center_list[[i]], f_center_list[[j]], -length(f_center_list[[i]]), length(f_center_list[[j]]), by=100)$dist_min
    }
  }
  return(list(clusters=clusters, f_center_var = f_center_var))
}





# Cluster curves (using GD) ----------------------------------------------------------

cluster_curves_gd = function(f_list, k, init_method='kmeans', step_size=0.05, MaxIter=100, seed=45, stopping_redu=0.01)
{
  if(init_method=='kmeans') init_index = init_kmeans(f_list, k, seed)
  # print(index)
  else if(init_method=='random') {set.seed(seed); init_index = sample(length(f_list), k)}

  
  f_center_list = f_list[init_index]
  
  n0_vec = rep(0, length(f_list))
  
  dist_curr = Inf
  dist_redu = Inf
  iter_count = 0
  while(dist_redu > stopping_redu && iter_count <= MaxIter)
  {
    iter_count = iter_count+1
    
    r = re_cluster_gd(f_list, f_center_list, n0_vec, step_size)
    clusters = r$clusters
    n0_vec = r$n0_vec
    
    r = re_center_gd(f_list, clusters, n0_vec, step_size)
    f_center_list = r$f_center_list
    dist_upd = r$mean_dist
    n0_vec = r$n0_vec
    
    dist_redu = (dist_curr-dist_upd) / (dist_upd)
    if (is.na(dist_redu)) dist_redu=0
    dist_curr = dist_upd
  }
  # clusters = re_cluster_gd(f_list, f_center_list, n0_vec, step_size)$clusters
  # print(iter_count)
  # print(dist_curr)
  # for (i in 1:length(f_center_list)) {
  #   plot(f_center_list[[i]],type='l')
  # }
  f_center_var = 0
  f_center_dist_min = Inf
  for (i in (1:(length(f_center_list)-1))) {
    for (j in ((i+1):length(f_center_list))) {
      dist = align_curves_gd(f_center_list[[i]], f_center_list[[j]], 0, step_size)$dist_min
      f_center_var = f_center_var+dist
      if (dist < f_center_dist_min) {
        f_center_dist_min = dist
      }
    }
  }
  return(list(clusters=clusters, f_center_var = f_center_var, f_center_dist_min = f_center_dist_min, f_center_list=f_center_list, n0_vec=n0_vec, init_index=init_index))
}


cluster_curves_gd_pp = function(f_list, k, init_method='kmeans', step_size=0.05, MaxIter=100, seed=45, stopping_redu=0.01)
{
  if(init_method=='kmeans') init_index = init_kmeans(f_list, k, seed, pp=TRUE)
  # print(index)
  else if(init_method=='random') {set.seed(seed); init_index = sample(length(f_list), k)}
  
  
  f_center_list = f_list[init_index]
  
  n0_vec = rep(0, length(f_list))
  
  dist_curr = Inf
  dist_redu = Inf
  iter_count = 0
  while(dist_redu > stopping_redu && iter_count <= MaxIter)
  {
    iter_count = iter_count+1
    
    r = re_cluster_gd(f_list, f_center_list, n0_vec, step_size, pp=TRUE)
    clusters = r$clusters
    n0_vec = r$n0_vec
    
    r = re_center_gd(f_list, clusters, n0_vec, step_size, pp=TRUE)
    f_center_list = r$f_center_list
    dist_upd = r$mean_dist
    n0_vec = r$n0_vec
    
    dist_redu = (dist_curr-dist_upd) / (dist_upd)
    if (is.na(dist_redu)) dist_redu=0
    dist_curr = dist_upd
  }
  # clusters = re_cluster_gd(f_list, f_center_list, n0_vec, step_size)$clusters
  # print(iter_count)
  # print(dist_curr)
  # for (i in 1:length(f_center_list)) {
  #   plot(f_center_list[[i]],type='l')
  # }
  f_center_var = 0
  f_center_dist_min = Inf
  if (length(f_center_list)>1){
    for (i in (1:(length(f_center_list)-1))) {
      for (j in ((i+1):length(f_center_list))) {
        dist = align_curves_gd(f_center_list[[i]], f_center_list[[j]], 0, step_size, pp=TRUE)$dist_min
        f_center_var = f_center_var+dist
        if (dist < f_center_dist_min) {
          f_center_dist_min = dist
        }
      }
    }
  }
  
  return(list(clusters=clusters, f_center_var = f_center_var, f_center_dist_min = f_center_dist_min, f_center_list=f_center_list, n0_vec=n0_vec, init_index=init_index))
}

# cluster_curves_gd(f_list, 2, 0.05)
# cluster_curves(f_list, 2)


# Test clustering result --------------------------------------------------

# {
#   x = seq(0,50,0.01)
#   f_1 = function(tau)
#     {return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,6+tau,2)))}
#   f_2 = function(tau)
#     {return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,9+tau,2)))}
#   f_3 = function(tau)
#     {return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,12+tau,2)))}
#   f_4 = function(tau)
#     {return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,3+tau,2)))}
#   
#   # f_center_1 = f_1(0)
#   # f_center_2 = f_2(0)
#   #
#   #
#   # N = length(f_center_1)
#   plot(f_1(0), type = 'l')
#   plot(f_2(0), type = 'l')
#   plot(f_3(0), type = 'l')
#   plot(f_4(0), type = 'l')
#   
#   
#   # What below is to generate shifted functions (shift to right)
#   clus_size_1 = 10
#   clus_size_2 = 10
#   clus_size_3 = 10
#   clus_size_4 = 10
#   set.seed(831); tau_vec_1 = runif(clus_size_1, min=0, max=20)
#   set.seed(93); tau_vec_2 = runif(clus_size_2, min=0, max=20)
#   set.seed(17); tau_vec_3 = runif(clus_size_3, min=0, max=20)
#   set.seed(65); tau_vec_4 = runif(clus_size_4, min=0, max=20)
#   
#   
#   f_list_1 = vector("list", clus_size_1)
#   f_list_2 = vector("list", clus_size_2)
#   f_list_3 = vector("list", clus_size_3)
#   f_list_4 = vector("list", clus_size_4)
#   
#   for (i in 1:clus_size_1) {
#     tau = tau_vec_1[i]
#     f_list_1[[i]] = f_1(tau)
#   }
#   for (i in 1:clus_size_2) {
#     tau = tau_vec_2[i]
#     f_list_2[[i]] = f_2(tau)
#   }
#   for (i in 1:clus_size_3) {
#     tau = tau_vec_3[i]
#     f_list_3[[i]] = f_3(tau)
#   }
#   for (i in 1:clus_size_4) {
#     tau = tau_vec_4[i]
#     f_list_4[[i]] = f_4(tau)
#   }
#   
#   f_list = c(f_list_1, f_list_2)
#   # f_list = c(f_list_1, f_list_2, f_list_3, f_list_4)
#   
#   # f_center_list = list(f_center_1, f_center_2)
#   
#   
#   ### check curves
#   # for (i in 1:clus_size_1) {
#   #   plot(f_list_1[[i]], type='l')
#   #   plot(f_list_2[[i]], type='l')
#   # }
#   
#   
#   cluster_curves(f_list, 4)
#   
# }










