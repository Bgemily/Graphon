source('~/Documents/Academic/SC/graphon/simulation/align_curves.R')


# Re-cluster --------------------------------------------------------------

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
  recur_count = 0
  while (dist_reduc>0.01 && recur_count <= MaxIter) {
    recur_count = recur_count + 1
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
  # print(recur_count)
  # plot(f_center, type = 'l')
  return(list(f_center=f_center, dist=dist_curr))
}


# clusters = list(c(1:5,16:20), c(6:15))
# re_center(f_list, clusters, mint_f_vec, maxt_f_vec)$mean_dist




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
  recur_count = 0
  while(dist_redu > stopping_redu && recur_count <= MaxIter)
  {
    recur_count = recur_count+1
    
    clusters = re_cluster(f_list, f_center_list, mint_f_vec, maxt_f_vec)
    r = re_center(f_list, clusters, mint_f_vec, maxt_f_vec)
    f_center_list = r$f_center_list
    dist_upd = r$mean_dist
    
    dist_redu = (dist_curr-dist_upd) / (dist_upd)
    if (is.na(dist_redu)) dist_redu=0
    dist_curr = dist_upd
  }
  clusters = re_cluster(f_list, f_center_list, mint_f_vec, maxt_f_vec)
  # print(recur_count)
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


# Test clustering result --------------------------------------------------

main = function()
{
  x = seq(0,50,0.01)
  f_1 = function(tau)
    return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,6+tau,2)))
  f_2 = function(tau)
    return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,9+tau,2)))
  f_3 = function(tau)
    return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,12+tau,2)))
  f_4 = function(tau)
    return(c(rep(0,length(x)),1/3*pnorm(x,5+tau,2)+2/3*pnorm(x,3+tau,2)))
  
  # f_center_1 = f_1(0)
  # f_center_2 = f_2(0)
  #
  #
  # N = length(f_center_1)
  plot(f_1(0), type = 'l')
  plot(f_2(0), type = 'l')
  plot(f_3(0), type = 'l')
  plot(f_4(0), type = 'l')
  
  
  # What below is to generate shifted functions (shift to right)
  clus_size_1 = 10
  clus_size_2 = 10
  clus_size_3 = 10
  clus_size_4 = 10
  set.seed(831); tau_vec_1 = runif(clus_size_1, min=0, max=20)
  set.seed(93); tau_vec_2 = runif(clus_size_2, min=0, max=20)
  set.seed(17); tau_vec_3 = runif(clus_size_3, min=0, max=20)
  set.seed(65); tau_vec_4 = runif(clus_size_4, min=0, max=20)
  
  
  f_list_1 = vector("list", clus_size_1)
  f_list_2 = vector("list", clus_size_2)
  f_list_3 = vector("list", clus_size_3)
  f_list_4 = vector("list", clus_size_4)
  
  for (i in 1:clus_size_1) {
    tau = tau_vec_1[i]
    f_list_1[[i]] = f_1(tau)
  }
  for (i in 1:clus_size_2) {
    tau = tau_vec_2[i]
    f_list_2[[i]] = f_2(tau)
  }
  for (i in 1:clus_size_3) {
    tau = tau_vec_3[i]
    f_list_3[[i]] = f_3(tau)
  }
  for (i in 1:clus_size_4) {
    tau = tau_vec_4[i]
    f_list_4[[i]] = f_4(tau)
  }
  
  # f_list = c(f_list_1, f_list_2)
  f_list = c(f_list_1, f_list_2, f_list_3, f_list_4)
  
  # f_center_list = list(f_center_1, f_center_2)
  
  
  ### check curves
  # for (i in 1:clus_size_1) {
  #   plot(f_list_1[[i]], type='l')
  #   plot(f_list_2[[i]], type='l')
  # }
  
  
  cluster_curves(f_list, 4)
  
}










