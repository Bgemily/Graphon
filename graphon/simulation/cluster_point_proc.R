source('~/Documents/Academic/SC/graphon/simulation/cluster_curves.R')


# Construct and cluster point process -------------------------------------

# main = function()
# {
  x = seq(0,50,0.01)
  
  num_nodes = 30
  clus_size_1 = 3
  clus_size_2 = 27
  # clus_size_3 = 10
  # clus_size_4 = 10
  set.seed(831); tau_vec_1 = runif(clus_size_1, min=0, max=20)
  set.seed(93); tau_vec_2 = runif(clus_size_2, min=0, max=20)
  # set.seed(17); tau_vec_3 = runif(clus_size_3, min=0, max=20)
  # set.seed(65); tau_vec_4 = runif(clus_size_4, min=0, max=20)
  
  N_list_1 = vector("list", clus_size_1)
  N_list_2 = vector("list", clus_size_2)
  set.seed(0); N_T_1_vec = sample(num_nodes, clus_size_1) # total number of edges
  set.seed(1); N_T_2_vec = sample(num_nodes, clus_size_2)
  f_list_1 = vector("list", clus_size_1)
  f_list_2 = vector("list", clus_size_2)
  
  for (i in 1:clus_size_1) {
    tau = tau_vec_1[i]
    N = N_T_1_vec[i]
    set.seed(4*i+1); N1 = rnorm(N%/%3, 5+tau, 2)
    set.seed(4*i+2); N2 = rnorm(N-N%/%3, 6+tau, 2)
    N_list_1[[i]] = c(N1,N2)
    
    f_empirical = function(t) sum(N_list_1[[i]]<=t)/length(N_list_1[[i]])
    f_list_1[[i]] = c(rep(0,length(x)), sapply(x, f_empirical))
  }
  for (i in 1:clus_size_2) {
    tau = tau_vec_2[i]
    N = N_T_2_vec[i]
    set.seed(4*i+3); N1 = rnorm(N%/%3, 5+tau, 2)
    set.seed(4*i+4); N2 = rnorm(N-N%/%3, 19+tau, 9)
    N_list_2[[i]] = c(N1,N2)
    
    f_empirical = function(t) sum(N_list_2[[i]]<=t)/length(N_list_2[[i]])
    f_list_2[[i]] = c(rep(0,length(x)), sapply(x, f_empirical))
  }
  f_list = c(f_list_1, f_list_2)
  cluster_curves(f_list, 2)
# }

  for (i in (1:9)) {
    plot(f_list[[i]], type = 'l', main = i)
  }

# Construct network -------------------------------------------------------
generate_network = function(SEED=0)
{
  grid_x = seq(0,2,by=.01)
  grid_y = seq(0,6,by=.01)
  grid = expand.grid(grid_x, grid_y)
  grid = as.matrix(grid)
  
  clus_size_1 = 3; clus_size_2 = 27
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2))
  centers = cbind(rep(1,clus_size_1), seq(0,6,length.out=clus_size_1))
  
  # clus_size_1 = 3; clus_size_2 = 3; clus_size_3 = 24
  # clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2), rep(3, clus_size_3))
  # centers = cbind(c(rep(1,clus_size_1), c(0.1,1.8,0.1)), c(seq(0,6,length.out=clus_size_1), seq(1,5,length.out=clus_size_2)))
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_2),])
  # nodes_mat = grid[sample(dim(grid)[1], clus_size_1+clus_size_2),]
  plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '')
  points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  
  # set.seed(42+SEED); 
  # nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_3),])
  # plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '')
  # points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  # points(nodes_mat[1:clus_size_1+clus_size_2,1], nodes_mat[1:clus_size_1+clus_size_2,2], col='blue')
  
  radius_thres = 1.5
  
  total_time = 50
  x = seq(0, total_time, 0.01)
  set.seed(98+SEED); tau_vec = runif(clus_size_1, 0, 30)
  # set.seed(98+SEED); tau_vec = c(runif(clus_size_1, 0, 10), runif(clus_size_2,10,20))
  
  # edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  # seed = 0+SEED
  # for (i in (1:length(clusters_vec))) {
  #   for (j in (1:length(clusters_vec))) {
  #     if (norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')>radius_thres)
  #       next
  #     if (clusters_vec[i]==clusters_vec[j]) {
  #       seed = seed+1; set.seed(seed)
  #       edge_time_mat[i,j] = runif(1, min=0, max=0.8*total_time)
  #     }
  #     else if (clusters_vec[i]==1 && clusters_vec[j]==3){
  #       seed = seed+1; set.seed(seed)
  #       tau = tau_vec[i]
  #       edge_time_mat[i,j] = rnorm(1, 5+tau, 1)
  #     }
  #     else if (clusters_vec[i]==2 && clusters_vec[j]==1)
  #   }
  # }
  
  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  for (i in (1:(length(clusters_vec)-1))) {
    for (j in ((i+1):length(clusters_vec))) {
      if (norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')>radius_thres)
        next
      if (clusters_vec[i]==clusters_vec[j]) {
        seed = seed+1; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.8*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
      else {
        seed = seed+1; set.seed(seed)
        if (clusters_vec[i]==1) tau = tau_vec[i]
        else tau = tau_vec[j]
        edge_time_mat[i,j] = rnorm(1, 5+tau, 1)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
    }
  }
  return(edge_time_mat)
}

get_emp_f_list = function(edge_time_mat)
{
  f_list = vector('list', dim(edge_time_mat)[1])
  for (i in (1:dim(edge_time_mat)[1])) {
    f_empirical = function(t) sum(edge_time_mat[i,]<=t)/sum(edge_time_mat[i,]!=Inf)
    f_list[[i]] = c(rep(0,length(x)), sapply(x, f_empirical))
  }
  return(f_list)
}


main = function(SEED, k=2)
{
  edge_time_mat = generate_network(SEED)
  f_list = get_emp_f_list(edge_time_mat)
  
  f_center_var = 0
  for (seed in (1:5+SEED)) {
    r = cluster_curves(f_list, k, seed=seed)
    if (r$f_center_var > f_center_var)
    {
      clusters = r$clusters
      f_center_var = r$f_center_var
    }
  }
  print(clusters)
  for (i in (1:length(f_list))) {
    plot(f_list[[i]], type='l', main=i)
  }
}


for (SEED in c(110)) {
  main(SEED)
}


# for (i in c(2,6,7,69,70)) {
#   plot(f_list[[i]], type = 'l', main = i)
# }
