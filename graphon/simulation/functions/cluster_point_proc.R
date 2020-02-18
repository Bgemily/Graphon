# source('~/Documents/Academic/SC/graphon/simulation/cluster_curves.R')
# source('./cluster_curves.R')


# Obtain empirical cdf from network -------------------------------------------------------
 

get_emp_f_list = function(edge_time_mat, t_vec)
{
  f_list = vector('list', dim(edge_time_mat)[1])
  for (i in (1:dim(edge_time_mat)[1])) {
    f_empirical = function(t) sum(edge_time_mat[i,]<=t)/sum(edge_time_mat[i,]<=max(t_vec))
    f_list[[i]] = c(rep(0,length(t_vec)), sapply(t_vec, f_empirical))
  }
  return(f_list)
}


# pdf
get_pp_f_list = function(edge_time_mat, t_vec, h=1)
{
  f_list = vector('list', dim(edge_time_mat)[1])
  for (i in (1:dim(edge_time_mat)[1])) {
    n = sum(edge_time_mat[i,]<=max(t_vec))
    f_density = function(t) 1/(n*h)*sum(sapply(edge_time_mat[i,], function(x)dnorm((t-x)/h, 0, 1)))
    f_list[[i]] = c(rep(0,length(t_vec)), sapply(t_vec, f_density))
}
  return(f_list)
}


# Network components ------------------------------------------------------

generate_grid = function(){
  grid_x = seq(0,1,by=.01)
  grid_y = seq(0,6,by=.01)
  grid = expand.grid(grid_x, grid_y)
  grid = as.matrix(grid)
  return(grid)
}


# Easiest case test --------------------------------------------------------------------

generate_network = function(SEED=0, total_time)
{
  t = seq(0, total_time, 0.05)
  
  grid = generate_grid()
  
  radius_thres = 1
  
  clus_size_1 = 4; clus_size_2 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2))
  set.seed(SEED+83); centers_x = runif(clus_size_1, 0.3,0.7)
  set.seed(SEED+724); centers_y = runif(clus_size_1, 0.8, 5.2)
  # centers = cbind(c(rep(0.5,clus_size_1)), c(seq(0.8,5.2,length.out=clus_size_1)))
  centers = cbind(centers_x, centers_y)
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_2),])
  
  set.seed(98+SEED); tau_vec = runif(clus_size_1, 0, 30)
  
  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  pdf_addend_list = vector("list", length(clusters_vec))
  for (i in (1:(length(clusters_vec)))) {
    for (j in ((1):length(clusters_vec))) {
      if (i==j || norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')>radius_thres)
        next
      if (clusters_vec[i]==2 && clusters_vec[j]==2 && edge_time_mat[i,j]==Inf) {
        seed = seed+1; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.8*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dunif(x, 0, 0.8*total_time))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dunif(x, 0, 0.8*total_time))))
      }
      else if(clusters_vec[i]==1 && clusters_vec[j]==2 && edge_time_mat[i,j]==Inf){
        seed = seed+1; set.seed(seed)
        tau = tau_vec[i]
        edge_time_mat[i,j] = tau + rnorm(1, 5, 1)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dnorm(x, tau+5, 1))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dnorm(x, tau+5, 1))))
      }
    }
  }
  
  pdf_list = list()
  for (i in 1:length(pdf_addend_list)) {
    if (dim(pdf_addend_list[[i]])[1]==1) pdf_list[[i]] = pdf_addend_list[[i]]
    pdf_list[[i]] = colMeans(pdf_addend_list[[i]])
  }
  
  
  return(list(edge_time_mat=edge_time_mat, nodes_mat=nodes_mat, tau_vec=tau_vec, pdf_list = pdf_list))
}



# Violate assumption but still okay ---------------------------------------


generate_network2 = function(SEED=0, total_time)
{
  t = seq(0, total_time, 0.05)
  
  grid = generate_grid()
  
  radius_thres1 = 2
  radius_thres2 = 1
  
  
  clus_size_1 = 4; clus_size_2 = 8; clus_size_3 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2), rep(3, clus_size_3))
  set.seed(SEED+487); centers_x = runif(clus_size_1+clus_size_2, 0.3, 0.7)
  set.seed(SEED+53); centers_y = runif(clus_size_1+clus_size_2, 0.8, 5.2)
  # centers = cbind(c(rep(0.5,clus_size_1), rep(c(0.7,0.3), clus_size_2/2)), c(seq(0.8,5.2,length.out=clus_size_1), seq(1,5,length.out=clus_size_2)))
  centers = cbind(centers_x, centers_y)
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_3),])
  
  set.seed(98+SEED); tau_vec = c(runif(clus_size_1, 40, 42), runif(clus_size_2,0,5))
  
  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  pdf_addend_list = vector("list", length(clusters_vec))
  for (i in (1:(length(clusters_vec)))) {
    for (j in ((1):length(clusters_vec))) {
      dij = norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')
      if (i==j || dij>radius_thres1)
        next
      if (clusters_vec[i]==3 && clusters_vec[j]==3 && dij<=radius_thres2 && edge_time_mat[i,j]==Inf) {
        seed = seed+124; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.6*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dunif(x, 0, 0.6*total_time))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dunif(x, 0, 0.6*total_time))))
      }
      else if ((clusters_vec[i]==2 && clusters_vec[j]==3) && dij<=radius_thres2 && edge_time_mat[i,j]==Inf) {
        seed = seed+165; set.seed(seed)
        tau = tau_vec[i]
        norm_mean = 5; norm_sd = 1
        edge_time_mat[i,j] = tau + rnorm(1, norm_mean, norm_sd)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dnorm(x, tau+norm_mean, norm_sd))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dnorm(x, tau+norm_mean, norm_sd))))
      }
      # else if (clusters_vec[i]==1 && clusters_vec[j]==2 && dij<=radius_thres1 && edge_time_mat[i,j]==Inf) {
      else if (clusters_vec[i]==1 && dij<=radius_thres1 && edge_time_mat[i,j]==Inf) {
        seed = seed+18; 
        tau = tau_vec[i]
        set.seed(seed)
        # edge_time_mat[i,j] = tau + runif(1, 0,6)
        unif_min = 0; unif_max = 6
        edge_time_mat[i,j] = tau + runif(1, unif_min, unif_max)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dunif(x, tau+unif_min, tau+unif_max))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dunif(x, tau+unif_min, tau+unif_max))))
      }
    }
  }
  pdf_list = list()
  for (i in 1:length(pdf_addend_list)) {
    if (dim(pdf_addend_list[[i]])[1]==1) pdf_list[[i]] = pdf_addend_list[[i]]
    pdf_list[[i]] = colMeans(pdf_addend_list[[i]])
  }
  return(list(edge_time_mat=edge_time_mat, nodes_mat=nodes_mat, tau_vec=tau_vec, pdf_list=pdf_list))
}




# Fail --------------------------------------------------------------------


generate_network3 = function(SEED=0, total_time)
{
  t = seq(0, total_time, 0.05)
  
  grid_x = seq(0,1,by=.01)
  grid_y = seq(0,6,by=.01)
  grid = expand.grid(grid_x, grid_y)
  grid = as.matrix(grid)
  
  radius_thres1 = 2
  radius_thres2 = 1
  
  clus_size_1 = 4; clus_size_2 = 8; clus_size_3 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2), rep(3, clus_size_3))
  set.seed(SEED+487); centers_x = runif(clus_size_1+clus_size_2, 0.3, 0.7)
  set.seed(SEED+53); centers_y = runif(clus_size_1+clus_size_2, 0.8, 5.2)
  # centers = cbind(c(rep(0.5,clus_size_1), rep(c(0.7,0.3), clus_size_2/2)), c(seq(0.8,5.2,length.out=clus_size_1), seq(1,5,length.out=clus_size_2)))
  centers = cbind(centers_x, centers_y)
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_3),])
  
  
  set.seed(98+SEED); tau_vec = c(runif(clus_size_1, 40, 42), runif(clus_size_2,0,30)) # this is the fail case, together with N(0,1)

  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  pdf_addend_list = vector("list", length(clusters_vec))
  for (i in (1:(length(clusters_vec)))) {
    for (j in ((1):length(clusters_vec))) {
      dij = norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')
      if (i==j || dij>radius_thres1)
        next
      if (clusters_vec[i]==3 && clusters_vec[j]==3 && dij<=radius_thres2 && edge_time_mat[i,j]==Inf) {
        seed = seed+124; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.6*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dunif(x, 0, 0.6*total_time))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dunif(x, 0, 0.6*total_time))))
      }
      else if ((clusters_vec[i]==2 && clusters_vec[j]==3) && dij<=radius_thres2 && edge_time_mat[i,j]==Inf) {
        seed = seed+165; set.seed(seed)
        tau = tau_vec[i]
        edge_time_mat[i,j] = tau + rnorm(1, 5,1)
        # edge_time_mat[i,j] = tau + rnorm(1, 5,4)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dnorm(x, tau+5, 1))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dnorm(x, tau+5, 1))))
      }
      # else if (clusters_vec[i]==1 && clusters_vec[j]==2 && dij<=radius_thres1 && edge_time_mat[i,j]==Inf) {
      else if (clusters_vec[i]==1 && dij<=radius_thres1 && edge_time_mat[i,j]==Inf) {
        seed = seed+18; 
        tau = tau_vec[i]
        
        set.seed(seed)
        edge_time_mat[i,j] = tau + runif(1, 0,6)
        edge_time_mat[j,i] = edge_time_mat[i,j]
        
        pdf_addend_list[[i]] = rbind(pdf_addend_list[[i]], (sapply(t, function(x)dunif(x, tau, tau+6))))
        pdf_addend_list[[j]] = rbind(pdf_addend_list[[j]], (sapply(t, function(x)dunif(x, tau, tau+6))))
      }
    }
  }
  pdf_list = list()
  for (i in 1:length(pdf_addend_list)) {
    if (dim(pdf_addend_list[[i]])[1]==1) pdf_list[[i]] = pdf_addend_list[[i]]
    pdf_list[[i]] = colMeans(pdf_addend_list[[i]])
  }
  return(list(edge_time_mat=edge_time_mat, nodes_mat=nodes_mat, tau_vec=tau_vec, pdf_list=pdf_list))
}




# Main --------------------------------------------------------------------

main = function(case, SEED, k=2, step_size=0.05)
{
  total_time = 50
  t = seq(0, total_time, 0.05)
  
  if (case==1) network = generate_network(SEED, total_time)
  if (case==2) network = generate_network2(SEED, total_time)
  if (case==3) network = generate_network3(SEED, total_time)
  
  edge_time_mat = network$edge_time_mat

  i=1
  while (i<=dim(edge_time_mat)[1]) {
    if(sum(edge_time_mat[i,]<Inf)==0)
    {
      cat('Deleting node', i, '\n')
      edge_time_mat = edge_time_mat[-i,-i]
    }
    else
      i = i+1
  }
  f_list = get_emp_f_list(edge_time_mat, t)
  
  f_center_var = 0
  f_center_dist_min = 0
  r_best = NULL
  for (seed in (1:3+SEED)) {
    # r = cluster_curves(f_list, k, seed=seed)
    r = cluster_curves_gd(f_list, k, seed=seed, step_size = step_size)
    # if (r$f_center_var > f_center_var)
    if (r$f_center_dist_min > f_center_dist_min)
    {
      clusters = r$clusters
      f_center_var = r$f_center_var
      f_center_dist_min = r$f_center_dist_min
      r_best = r
    }
    # print(r$f_center_var)
    # print(r$f_center_dist_min)
  }
  # print(clusters)
  return(list(network=network, f_list=f_list, f_center_list=r_best$f_center_list, clusters=r_best$clusters, n0_ve=r_best$n0_vec, f_center_var=r_best$f_center_var, init_index=r_best$init_index))
}



main_pp = function(case, SEED, k=2, step_size=0.05, h=1)
{
  total_time = 50
  t = seq(0, total_time, 0.05)
  
  if (case==1) network = generate_network(SEED, total_time)
  if (case==2) network = generate_network2(SEED, total_time)
  if (case==3) network = generate_network3(SEED, total_time)
  
  edge_time_mat = network$edge_time_mat
  
  i=1
  while (i<=dim(edge_time_mat)[1]) {
    if(sum(edge_time_mat[i,]<Inf)==0)
    {
      cat('Deleting node', i, '\n')
      edge_time_mat = edge_time_mat[-i,-i]
    }
    else
      i = i+1
  }
  
  
  # f_list = get_emp_f_list(edge_time_mat, t)
  f_list = get_pp_f_list(edge_time_mat, t, h=h)
  
  f_center_var = 0
  f_center_dist_min = 0
  r_best = NULL
  for (seed in (1:3+SEED)) {
    r = cluster_curves_gd_pp(f_list, k, seed=seed, step_size = step_size)
    # if (r$f_center_var > f_center_var)
    if (r$f_center_dist_min > f_center_dist_min)
    {
      clusters = r$clusters
      f_center_var = r$f_center_var
      f_center_dist_min = r$f_center_dist_min
      r_best = r
    }
    # print(r$f_center_var)
    # print(r$f_center_dist_min)
  }
  # print(clusters)
  return(list(network=network, f_list=f_list, f_center_list=r_best$f_center_list, clusters=r_best$clusters, n0_ve=r_best$n0_vec, f_center_var=r_best$f_center_var, init_index=r_best$init_index))
}






