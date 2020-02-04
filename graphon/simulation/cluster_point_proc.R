source('~/Documents/Academic/SC/graphon/simulation/cluster_curves.R')


# Obtain empirical cdf from network -------------------------------------------------------
 

get_emp_f_list = function(edge_time_mat, t)
{
  f_list = vector('list', dim(edge_time_mat)[1])
  for (i in (1:dim(edge_time_mat)[1])) {
    f_empirical = function(t) sum(edge_time_mat[i,]<=t)/sum(edge_time_mat[i,]!=Inf)
    f_list[[i]] = c(rep(0,length(t)), sapply(t, f_empirical))
  }
  return(f_list)
}


# Easiest case test --------------------------------------------------------------------

generate_network = function(SEED=0, total_time)
{
  grid_x = seq(0,1,by=.01)
  grid_y = seq(0,6,by=.01)
  grid = expand.grid(grid_x, grid_y)
  grid = as.matrix(grid)
  
  radius_thres = 1
  
  clus_size_1 = 4; clus_size_2 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2))
  centers = cbind(c(rep(0.5,clus_size_1)), c(seq(0.8,5.2,length.out=clus_size_1)))
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_2),])
  # dev.new(width=2,height=6,noRStudioGD = T)
  plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '', xlim=c(0,1))
  points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  
  # plot the circle
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points(x_center+radius_thres*cos(angel), y_center+radius_thres*sin(angel), cex=0.1, col='red')
  
  
  set.seed(98+SEED); tau_vec = runif(clus_size_1, 0, 30)
  
  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  for (i in (1:(length(clusters_vec)))) {
    for (j in ((1):length(clusters_vec))) {
      if (i==j || norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')>radius_thres)
        next
      if (clusters_vec[i]==2 && clusters_vec[j]==2) {
        seed = seed+1; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.8*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
      else if(clusters_vec[i]==1 && clusters_vec[j]==2){
        seed = seed+1; set.seed(seed)
        tau = tau_vec[i]
        edge_time_mat[i,j] = tau + rnorm(1, 5, 1)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
    }
  }
  return(list(edge_time_mat=edge_time_mat, nodes_mat=nodes_mat))
}


main = function(SEED, k=2, step_size=0.05)
{
  total_time = 50
  t = seq(0, total_time, 0.01)
  
  # edge_time_mat = generate_network2(SEED, total_time)
  edge_time_mat = generate_network(SEED, total_time)$edge_time_mat
  f_list = get_emp_f_list(edge_time_mat, t)
  
  f_center_var = 0
  r_best = NULL
  for (seed in (1:3+SEED)) {
    r = cluster_curves_gd(f_list, k, seed=seed, step_size = step_size)
    if (r$f_center_var > f_center_var)
    {
      clusters = r$clusters
      f_center_var = r$f_center_var
      r_best = r
    }
    print(r$f_center_var)
  }
  print(clusters)
  return(list(f_list=f_list, f_center_list=r_best$f_center_list, clusters=r_best$clusters, n0_ve=r_best$n0_vec, f_center_var=r_best$f_center_var))
}


SEED_vec = seq(108,117,length.out=10)
results = vector("list", 0)
for (SEED in SEED_vec) {
  results[[as.character(SEED)]]=main(SEED, k=2, step_size = 0.01)
}




# Violate assumption but still okay ---------------------------------------


generate_network2 = function(SEED=0, total_time)
{
  grid_x = seq(0,1,by=.01)
  grid_y = seq(0,6,by=.01)
  grid = expand.grid(grid_x, grid_y)
  grid = as.matrix(grid)
  
  radius_thres1 = 2
  radius_thres2 = 1
  
  
  clus_size_1 = 4; clus_size_2 = 8; clus_size_3 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2), rep(3, clus_size_3))
  centers = cbind(c(rep(0.5,clus_size_1), rep(c(0.7,0.3), clus_size_2/2)), c(seq(0.8,5.2,length.out=clus_size_1), seq(1,5,length.out=clus_size_2)))
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_3),])
  # dev.new(width=2,height=6,noRStudioGD = T)
  plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '', xlim=c(0,1))
  points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  points(nodes_mat[1:clus_size_2+clus_size_1,1], nodes_mat[1:clus_size_2+clus_size_1,2], col='blue')
  
  # plot the circles
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points(x_center+radius_thres1*cos(angel), y_center+radius_thres1*sin(angel), cex=0.1, col='red')
  x_center = centers[1+clus_size_1,1]; y_center = centers[1+clus_size_1,2]
  points(x_center+radius_thres2*cos(angel), y_center+radius_thres2*sin(angel), cex=0.1, col='blue')
  
  
  set.seed(98+SEED); tau_vec = c(runif(clus_size_1, 40, 42), runif(clus_size_2,0,10))
  
  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  for (i in (1:(length(clusters_vec)))) {
    for (j in ((1):length(clusters_vec))) {
      dij = norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')
      if (i==j || dij>radius_thres1)
        next
      if (clusters_vec[i]==3 && clusters_vec[j]==3 && dij<=radius_thres2) {
        seed = seed+124; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.6*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
      else if ((clusters_vec[i]==2 && clusters_vec[j]==3) && dij<=radius_thres2) {
        seed = seed+165; set.seed(seed)
        tau = tau_vec[i]
        edge_time_mat[i,j] = tau + rnorm(1, 5,4)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
      else if (clusters_vec[i]==1 && clusters_vec[j]==2 && dij<=radius_thres1) {
        seed = seed+18; 
        tau = tau_vec[i]
        set.seed(seed)
        edge_time_mat[i,j] = tau + runif(1, 0,6)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
    }
  }
  return(list(edge_time_mat=edge_time_mat, nodes_mat=nodes_mat))
}



main2 = function(SEED, k=2, step_size=0.05)
{
  total_time = 50
  t = seq(0, total_time, 0.01)
  
  edge_time_mat = generate_network2(SEED, total_time)$edge_time_mat
  f_list = get_emp_f_list(edge_time_mat, t)
  
  f_center_var = 0
  r_best = NULL
  for (seed in (1:3+SEED)) {
    # r = cluster_curves(f_list, k, seed=seed)
    r = cluster_curves_gd(f_list, k, seed=seed, step_size = step_size)
    if (r$f_center_var > f_center_var)
    {
      clusters = r$clusters
      f_center_var = r$f_center_var
      r_best = r
    }
    print(r$f_center_var)
  }
  print(clusters)
  return(list(f_list=f_list, f_center_list=r_best$f_center_list, clusters=r_best$clusters, n0_ve=r_best$n0_vec, f_center_var=r_best$f_center_var))
}


# SEED_vec = c(109)
SEED_vec = seq(108,117,length.out=10)
results2 = vector("list", 0)
for (SEED in SEED_vec) {
  results2[[as.character(SEED)]]=main2(SEED, k=3, step_size = 0.01)
}






total_time = 50
t = seq(0, total_time, 0.01)
edge_time_mat = generate_network2(117, total_time)$edge_time_mat
f_list = get_emp_f_list(edge_time_mat, t)
for (i in c(13,31,35,44)) {
  plot(f_list[[i]], type='l', main=i)
}


# Fail --------------------------------------------------------------------


generate_network3 = function(SEED=0, total_time)
{
  grid_x = seq(0,1,by=.01)
  grid_y = seq(0,6,by=.01)
  grid = expand.grid(grid_x, grid_y)
  grid = as.matrix(grid)
  
  radius_thres1 = 2
  radius_thres2 = 1
  
  
  clus_size_1 = 4; clus_size_2 = 8; clus_size_3 = 44
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2), rep(3, clus_size_3))
  centers = cbind(c(rep(0.5,clus_size_1), rep(c(0.7,0.3), clus_size_2/2)), c(seq(0.8,5.2,length.out=clus_size_1), seq(1,5,length.out=clus_size_2)))
  
  set.seed(42+SEED);
  nodes_mat = rbind(centers, grid[sample(dim(grid)[1], clus_size_3),])
  # par(pin=c(1,6)*.4)
  plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '', xlim=c(0,1))
  points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  points(nodes_mat[1:clus_size_2+clus_size_1,1], nodes_mat[1:clus_size_2+clus_size_1,2], col='blue')
  
  # plot the circle
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points(x_center+radius_thres1*cos(angel), y_center+radius_thres1*sin(angel), cex=0.1, col='red')
  x_center = centers[1+clus_size_1,1]; y_center = centers[1+clus_size_1,2]
  points(x_center+radius_thres2*cos(angel), y_center+radius_thres2*sin(angel), cex=0.1, col='blue')
  
  
  set.seed(98+SEED); tau_vec = c(runif(clus_size_1, 40, 42), runif(clus_size_2,0,10))
  
  edge_time_mat = matrix(Inf, nrow=length(clusters_vec), ncol=length(clusters_vec))
  seed = 0+SEED
  for (i in (1:(length(clusters_vec)))) {
    for (j in ((1):length(clusters_vec))) {
      dij = norm(t(nodes_mat[i,]-nodes_mat[j,]), 'f')
      if (i==j || dij>radius_thres1)
        next
      if (clusters_vec[i]==3 && clusters_vec[j]==3 && dij<=radius_thres2) {
        seed = seed+124; set.seed(seed)
        edge_time_mat[i,j] = runif(1, min=0, max=0.6*total_time)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
      else if ((clusters_vec[i]==2 && clusters_vec[j]==3) && dij<=radius_thres2) {
        seed = seed+165; set.seed(seed)
        tau = tau_vec[i]
        edge_time_mat[i,j] = tau + rnorm(1, 5,4)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
      # else if (clusters_vec[i]==1 && dij<=radius_thres2) {
      else if (clusters_vec[i]==1 && clusters_vec[j]==2 && dij<=radius_thres1) {
        seed = seed+18; 
        tau = tau_vec[i]
        # set.seed(seed+.9); edge_prob = runif(1)
        # if(edge_prob<0.7)
        # {
        #   set.seed(seed)
        #   edge_time_mat[i,j] = tau + rnorm(1, 5,2)
        #   edge_time_mat[j,i] = edge_time_mat[i,j]
        # }
        # else
        # {
        #   set.seed(seed)
        #   edge_time_mat[i,j] = tau-35 + rnorm(1, 5,2)
        #   edge_time_mat[j,i] = edge_time_mat[i,j]
        # }
        set.seed(seed)
        edge_time_mat[i,j] = tau + runif(1, 0,6)
        edge_time_mat[j,i] = edge_time_mat[i,j]
      }
    }
  }
  return(list(edge_time_mat=edge_time_mat, nodes_mat=nodes_mat))
}



main3 = function(SEED, k=2, step_size=0.05)
{
  total_time = 50
  t = seq(0, total_time, 0.01)
  
  edge_time_mat = generate_network3(SEED, total_time)$edge_time_mat
  f_list = get_emp_f_list(edge_time_mat, t)
  
  # f_list = f_list[1:12]
  f_center_var = 0
  for (seed in (1:3+SEED)) {
    # r = cluster_curves(f_list, k, seed=seed)
    r = cluster_curves_gd(f_list, k, seed=seed, step_size = step_size)
    if (r$f_center_var > f_center_var)
    {
      clusters = r$clusters
      f_center_var = r$f_center_var
    }
    print(r$f_center_var)
  }
  print(clusters)
  # for (i in (1:length(f_list))) {
  #   plot(f_list[[i]], type='l', main=i)
  # }
}


# main3(SEED, k=3, step_size = 0.01)

for (SEED in seq(108,117,length.out=10)) {
  main3(SEED, k=3, step_size = 0.01)
}



total_time = 50
t = seq(0, total_time, 0.01)

edge_time_mat = generate_network3(117, total_time)$edge_time_mat
f_list = get_emp_f_list(edge_time_mat, t)
for (i in c(1:20)) {
  plot(f_list[[i]], type='l', main=i)
}






# Visualization -----------------------------------------------------------

total_time = 50
t = seq(0, total_time, 0.01)

r = results$`108`
edge_time_mat = generate_network(108, total_time)$edge_time_mat
nodes_mat = generate_network(108, total_time)$nodes_mat


### plot the growing network
{
    library(colorRamps)
    library(grDevices)
    library(fields)
    # colorbar = colorRamp(c(rgb(1,0,0,0), rgb(1,0,0,1)), alpha=T)
    colorbar = cm.colors(51)
    colorbar = blue2red(50)
    
    
    dev.new(width=3, height=6, noRStudioGD = T)
    
    par(oma=c( 0,0,0,4))
    plot(nodes_mat[,1], nodes_mat[,2], cex = .5, xlab='', ylab = '', xlim=c(0,1), ylim=c(0,6))
    edge_time_round_mat = round(edge_time_mat)
    for (t in 1:50) {
      node_index_mat = which(edge_time_round_mat==t, arr.ind=T)
      if (dim(node_index_mat)[1]==0) {
        next
      }
      for (i in 1:dim(node_index_mat)[1]) {
        line = node_index_mat[i,]
        lines(nodes_mat[line,1], nodes_mat[line,2], col=colorbar[t+1])
      }
    }
    points(nodes_mat[1:4,1], nodes_mat[1:4,2], col='red')
    
    par(oma=c( 0,0,0,0))
    image.plot(legend.only = T,legend.lab = "time", col=colorbar, zlim=c(0,50))

}

### plot original cdf's
f_list = r$f_list

colors = c(rep('red', 4), rep('black', 46))
dev.new(width=6, height=4, noRStudioGD = T)
par(mfrow=c(1,2))
plot(t, tail(f_list[[1]], length(t)), type='l',col=colors[1], xlab='time', ylab='f(t)', main='Cluster 1')
for (i in 2:4) {
  lines(t, tail(f_list[[i]], length(t)), col=colors[i])
}
plot(t, tail(f_list[[5]], length(t)), type='l',col=colors[5], xlab='time', ylab='f(t)', main='Cluster 2')
for (i in 6:50) {
  lines(t, tail(f_list[[i]], length(t)), col=colors[i])
}


# colors = c(rep('red', 4), rep('blue', 8), rep('black', 46))
# dev.new(width=8, height=4, noRStudioGD = T)
# par(mfrow=c(1,3))
# plot(t, tail(f_list[[1]], length(t)), type='l',col=colors[1], xlab='time', ylab='f(t)', main='Cluster 1')
# for (i in 2:4) {
#   lines(t, tail(f_list[[i]], length(t)), col=colors[i])
# }
# plot(t, tail(f_list[[5]], length(t)), type='l',col=colors[5], xlab='time', ylab='f(t)', main='Cluster 2')
# for (i in 6:12) {
#   lines(t, tail(f_list[[i]], length(t)), col=colors[i])
# }
# plot(t, tail(f_list[[13]], length(t)), type='l',col=colors[13], xlab='time', ylab='f(t)', main='Cluster 3')
# for (i in 14:58) {
#   lines(t, tail(f_list[[i]], length(t)), col=colors[i])
# }


### plot aligned cdf's as well as the centers
n0_vec = r$n0_ve
f_center_list = r$f_center_list
clusters = r$clusters

colors = c(rep('red', 4), rep('black', 46))
dev.new(width=6, height=4, noRStudioGD = T)
par(mfrow=c(1,2))
plot(t, tail(shift(f_center_list[[1]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='black', xlab='time', ylab='f(t)', main='')
for (i in clusters[[1]]) {
  lines(t, tail(shift(f_list[[i]], n0_vec[i]), length(t)), col=colors[i], lwd = 0.5)
}
plot(t, tail(shift(f_center_list[[2]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='red', xlab='time', ylab='f(t)', main='')
for (i in clusters[[2]]) {
  lines(t, tail(shift(f_list[[i]], n0_vec[i]), length(t)), col=colors[i], lwd = 0.5)
}


# colors = c(rep('red', 4), rep('blue', 8), rep('black', 46))
# dev.new(width=8, height=4, noRStudioGD = T)
# par(mfrow=c(1,3))
# plot(t, tail(shift(f_center_list[[1]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='red', xlab='time', ylab='f(t)', main='Cluster 1')
# for (i in clusters[[1]]) {
#   lines(t, tail(shift(f_list[[i]], n0_vec[i]), length(t)), col=colors[i], lwd = 0.5)
# }
# plot(t, tail(shift(f_center_list[[2]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='black', xlab='time', ylab='f(t)', main='Cluster 2')
# for (i in clusters[[2]]) {
#   lines(t, tail(shift(f_list[[i]], n0_vec[i]), length(t)), col=colors[i], lwd = 0.5)
# }
# plot(t, tail(shift(f_center_list[[3]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='blue', xlab='time', ylab='f(t)', main='Cluster 3')
# for (i in clusters[[3]]) {
#   lines(t, tail(shift(f_list[[i]], n0_vec[i]), length(t)), col=colors[i], lwd = 0.5)
# }


### plot centers
dev.new(width=6, height=6, noRStudioGD = T)
plot(t, tail(shift(f_center_list[[1]], 0), length(t)), ylim = c(0,1), type='l', lty='dashed', lwd=1.5, col='black', xlab='time', ylab='f(t)')
lines(t, tail(shift(f_center_list[[2]], 0), length(t)), type='l', lty='solid', lwd=1.5, col='red')
legend("topleft", legend=c("Cluster 1", "Cluster 2"), col=c("black", "red"), lty=c('dashed','dashed'), cex=0.8)



# dev.new(width=6, height=6, noRStudioGD = T)
# plot(t, tail(shift(f_center_list[[1]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='red', xlab='time', ylab='f(t)')
# lines(t, tail(shift(f_center_list[[2]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='black')
# lines(t, tail(shift(f_center_list[[3]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='blue')
# legend("topleft", legend=c("Cluster 1", "Cluster 2", 'Cluster 3'), col=c("red", "black", "blue"), lty=c('dashed','dashed','dashed'), cex=0.8)
# 





