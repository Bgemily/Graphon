

# case 1 --------------------------------------------------------------------

generate_network1 = function(SEED=0, total_time=50)
{
  
  set.seed(SEED)
  
  t_vec = seq(0, total_time, 0.05)
  dist_thres = 2
  N_clus = 2
  clus_size_vec = c(50,50)
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)
  
  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_vec = c(runif(clus_size_vec[1], 0, 20), rep(0, clus_size_vec[2]))
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  # pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=1, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  # pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  # pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1));  
  # pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  # pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
      clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
      if (q==l){
        tmp = mapply(min, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
        tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
      }
      else if (q<l)
        tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
      
      else
        tau_mat[clus_q, clus_l] = t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q))
    }
  }
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
  
}



# case 2 ---------------------------------------


generate_network2 = function(SEED=0, total_time=50)
{
  set.seed(SEED)
  
  t_vec = seq(0, total_time, 0.05)
  dist_thres = 2
  N_clus = 3
  clus_size_vec = c(30,30,30)
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)
  
  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_vec = c(runif(clus_size_vec[1], 25, 30), runif(clus_size_vec[2],0,10), rep(0, clus_size_vec[3]))
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=10, rate=2), random = function(n) rgamma(n,10,2)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=10, rate=1), random = function(n) rgamma(n,10,1)); 
  pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dgamma(x, shape=20, rate=4), random = function(n) rgamma(n,20,4));  
  pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
      clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
      if (q==l){
        tmp = mapply(min, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
        # tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
        tau_mat[clus_q, clus_l] = 0
      }
      else if (q<l)
        tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
      
      else
        tau_mat[clus_q, clus_l] = t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q))
    }
  }
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
  
}



# case 3 ------------------------------------------------------------------


generate_network3 = function(SEED=0, total_time=50)
{
  set.seed(SEED)
  
  t_vec = seq(0, total_time, 0.05)
  dist_thres = 2
  N_clus = 3
  clus_size_vec = c(30,30,30)
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)
  
  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  tau_mat[clus_true[[1]], clus_true[[1]]] = 0
  tau_mat[clus_true[[2]], clus_true[[2]]] = 0
  tau_mat[clus_true[[3]], clus_true[[3]]] = 0
  tau_mat[clus_true[[1]], clus_true[[2]]] = runif(clus_size_vec[1]*clus_size_vec[2], 40, 42)
  tau_mat[clus_true[[1]], clus_true[[3]]] = runif(clus_size_vec[1]*clus_size_vec[3], 30, 32)
  tau_mat[clus_true[[2]], clus_true[[3]]] = runif(clus_size_vec[2]*clus_size_vec[3], 0, 10)
  tau_mat[lower.tri(tau_mat)] = t(tau_mat)[lower.tri(tau_mat)]
  
  if(!isSymmetric(tau_mat))
    stop("tau_mat should be symmetric.")
  
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=2, rate=1), random = function(n) rgamma(n,2,1)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dgamma(x, shape=10, rate=2), random = function(n) rgamma(n,10,2));  
  pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
  
}


# case 4 --------------------------------------------------------------------


generate_network4 = function(SEED=0, total_time=50)
{
  set.seed(SEED)
  
  t_vec = seq(0, total_time, 0.05)
  dist_thres = 2
  N_clus = 3
  clus_size_vec = c(50,50,50)
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)

  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_vec = c(runif(clus_size_vec[1], 40, 42), runif(clus_size_vec[2],0,2), rep(0, clus_size_vec[3]))
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1));  
  pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1))
  pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
      clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
      if (q==l){
        tmp = mapply(min, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
        tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
      }
      else if (q<l)
        tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
      
      else
        tau_mat[clus_q, clus_l] = t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q))
    }
  }
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
}



