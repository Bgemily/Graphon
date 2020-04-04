
rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Setup  ------------------------------------------------------------------

total_time = 50
time_unit = 0.05
t = seq(0, total_time, 0.05)
pp = TRUE

r = results2[[7]]

f_list = r$f_list
f_center_list = r$f_center_list
clusters = r$clusters
n0_mat = r$n0_mat
dist_mat = r$dist_mat
n0_vec = r$n0_vec
network = r$network
edge_time_mat = network$edge_time_mat
nodes_mat = network$nodes_mat
pdf_list = network$pdf_list
membership_true = network$membership_true
tau_vec = network$tau_vec


membership_est = unlist(clusters)
for (j in 1:length(clusters)) {
  membership_est[clusters[[j]]] = j
}

par(mar=c(2.5,2.5,.5,.5))

# ARI ---------------------------------------------------------------------

jitter_boxplot = function(ARI, group=1)
{  
  library(ggplot2)
  data_frame = data.frame(ARI, group=group)
  ggplot(data_frame, aes(x=group, y=ARI,group=group)) + 
    geom_boxplot(outlier.shape=NA)+
    # geom_violin()
    geom_jitter(position=position_jitter(width=.1, height=0))
}

jitter_boxplot(ARI=c(ARI3_old, ARI3),group=c(rep('',length(ARI3_old)),rep('cont.cluster',length(ARI3))))
jitter_boxplot(ARI3_cont)

membership_true1 = results1[[1]]$network$membership_true
ARI1 = get_ARI(membership_true1, results1, length(membership_true1))


membership_true2 = results2[[1]]$network$membership_true
ARI2 = get_ARI(membership_true2, results2, length(membership_true2))
ARI2_old = get_ARI(membership_true2, results2, length(membership_true2), clusters_old = TRUE)



membership_true3 = results3[[1]]$network$membership_true
ARI3 = get_ARI(membership_true3, results3, length(membership_true3))
ARI3_old = get_ARI(membership_true3, results3, length(membership_true3), clusters_old = TRUE)
ARI3_cont = get_ARI(membership_true3, results3_cont, length(membership_true3))



membership_true4 = results4[[1]]$network$membership_true
ARI4 = get_ARI(membership_true4, results4, length(membership_true4))


# plot connecting pattern matrix ------------------------------------------
{
true_pattern_matrix = vector("list",2)
true_pattern_matrix[[1]] = vector("list",3)
true_pattern_matrix[[2]] = vector("list",3)
# true_pattern_matrix[[3]] = vector("list",3)
true_pattern_matrix[[1]][[1]] = function(x) {mean=5; dnorm(x, mean, 1)}
true_pattern_matrix[[1]][[2]] = function(x) {mean=5; dnorm(x, mean, 1)}
# true_pattern_matrix[[1]][[3]] = function(x) {tau=40; dunif(x, tau, tau+6)}
true_pattern_matrix[[2]][[2]] = function(x) {tau=40; dunif(x, 0, 50)}
# true_pattern_matrix[[2]][[3]] = function(x) {mean=5; dnorm(x, mean, 0.5)}
# true_pattern_matrix[[3]][[3]] = function(x) {tau=0; dunif(x, tau, tau+30)}
true_pattern_matrix[[2]][[1]] = true_pattern_matrix[[1]][[2]]
# true_pattern_matrix[[3]][[1]] = true_pattern_matrix[[1]][[3]]
# true_pattern_matrix[[3]][[2]] = true_pattern_matrix[[2]][[3]]
}


plot_pattern_matrix = function(clusters_overclus, edge_time_mat, n0_mat)
{
  k = length(clusters_overclus)
  n0_vec_overclus = rep(0, nrow(edge_time_mat))
  for (i in 1:k) {
    n0_vec_overclus[clusters_overclus[[i]]] = n0_mat[clusters_overclus[[i]], clusters_overclus[[i]][1]] - min(n0_mat[clusters_overclus[[i]], clusters_overclus[[i]][1]])
  }
  
  pdf_array = get_pdf_array_(clusters_overclus, edge_time_mat, n0_vec_overclus)
  
  par(mfrow = c(k,k))
  for (i in 1:k) {
    for (j in 1:k) {
      # edge_time_submat = edge_time_mat[clusters_overclus[[i]], clusters_overclus[[j]]]
      # 
      # tau_i_vec = time_unit * n0_vec_overclus[clusters_overclus[[i]]]
      # tau_j_vec = time_unit * n0_vec_overclus[clusters_overclus[[j]]]
      # # tau_ij_mat = 1/2*(matrix(tau_i_vec, length(tau_i_vec), length(tau_j_vec)) + t(matrix(tau_j_vec, length(tau_j_vec), length(tau_i_vec))))
      # tau_ij_mat = matrix(tau_i_vec, length(tau_i_vec), length(tau_j_vec))
      # 
      # edge_time_submat = matrix(edge_time_submat-tau_ij_mat, nrow=1)
      # if (sum(edge_time_submat<Inf)==0) pdf = rep(0, 2*length(t))
      # else pdf = get_pp_f_list(edge_time_submat, t, h=1)[[1]]
      
      pdf = pdf_array[i,j,]
      plot(t, tail(pdf, length(t)), type='l', ylim = c(0,0.35))
      
      # lines(t, sapply(t,true_pattern_matrix[[i]][[j]]), col=2)
    }
  }
  
  par(mfrow = c(1,1)); par(mfrow = c(1,1))
  return()
}


plot(cmdscale(dist_mat), pch=membership_true, col=membership_est)
plot_pattern_matrix(clusters, edge_time_mat, n0_mat)


k=5
W = exp(-dist_mat^2/median(dist_mat)^2)
membership_overclus = spectral_clustering(W, k)
plot(cmdscale(dist_mat), pch=membership_true, col=membership_overclus)

clusters_overclus = list()
for (i in 1:k) {
  clusters_overclus[[i]] = which(membership_overclus==i)
}

clusters_merged = merge_clusters(clusters_overclus, edge_time_mat, n0_mat, 3)

membership_merged = unlist(clusters_merged)
for (j in 1:length(clusters_merged)) {
  membership_merged[clusters_merged[[j]]] = j
}
plot(cmdscale(dist_mat), pch=membership_true, col=membership_merged)


plot_pattern_matrix(clusters_overclus, edge_time_mat, n0_mat)
plot_pattern_matrix(clusters_merged, edge_time_mat, n0_mat)


# overcluster and merge ---------------------------------------------------

k_trueclus = 3
k_overclus = 6
results = results4[1:10]
ARI = ARI4[1:10]

membership_true = results[[1]]$network$membership_true


results_merge = vector("list", length(results)) 
for (i in 1:length(results)) {
  dist_mat = results[[i]]$dist_mat
  edge_time_mat = results[[i]]$network$edge_time_mat
  n0_mat = results[[i]]$n0_mat
  
  W = exp(-dist_mat^2/median(dist_mat)^2)
  membership_overclus = spectral_clustering(W, k_overclus)
  
  clusters_overclus = list()
  for (j in 1:k_overclus) {
    clusters_overclus[[j]] = which(membership_overclus==j)
  }
  results_merge[[i]] = list(clusters=merge_clusters(clusters_overclus, edge_time_mat, n0_mat, k_trueclus))
}

ARI_merge = get_ARI(membership_true, results_merge)
jitter_boxplot(c(ARI, ARI_merge), group = c(rep(paste('k=',k_trueclus,sep=''),length(ARI)),rep(paste('k=',k_overclus,' and merge',sep=''),length(ARI_merge))))

# plot(ARI, ARI_merge,xlim=c(0,1),ylim=c(0,1))
# abline(a=0,b=1,col='red')



# plot pairwise distance by multidimensional scaling --------------------------------------------------


# cluster result
plot(cmdscale(dist_mat), pch=membership_true, col=membership_est)

W = exp(-dist_mat^2/median(dist_mat)^2)
image(W)

# distance between true pdf's and cluster result
dist_mat_pp = sqrt(pairwise_dist_mat(pdf_list, pp=TRUE, n0_mat=r$n0_mat)$dist_mat)
plot(cmdscale(dist_mat_pp), pch=membership_true, col=membership_est)



# Plot nodes locations with cluster results -------------------------------


membership_est = unlist(clusters)
for (j in 1:length(clusters)) {
  membership_est[clusters[[j]]] = j
}


{ radius_thres1 = 2
  radius_thres2 = 2
  
  clus_size_1 = 30; clus_size_2 = 30; clus_size_3 = 30
  centers = nodes_mat[1:(clus_size_1+clus_size_2),]
  
  dev.new(width=6,height=1.5,noRStudioGD = T)
  par(mar = c(2.5,2.5,1,1))
  plot( nodes_mat[,2], nodes_mat[,1], pch = membership_true, col = membership_est, cex = 1, xlab='', ylab = '', xlim=c(0,6), ylim=c(0,1))
  
  
  # # plot the neighboor area
  # angel = seq(0, 2*pi, length.out=200)
  # x_center = centers[1,1]; y_center = centers[1,2]
  # points(y_center, x_center, pch = 21, bg='darkgray')
  # points( y_center+radius_thres1*sin(angel), x_center+radius_thres1*cos(angel), cex=0.1, col='darkgray')
  # 
  # x_center = centers[1+clus_size_1,1]; y_center = centers[1+clus_size_1,2]
  # points(y_center, x_center, pch = 24, bg='darkgray')
  # points(y_center+radius_thres2*sin(angel), x_center+radius_thres2*cos(angel), cex=0.1, col='darkgray')
}

{
  radius_thres1 = 1

  clus_size_1 = 4; clus_size_2 =  46
  centers = nodes_mat[1:(clus_size_1),]
  
  dev.new(width=6,height=1.5,noRStudioGD = T)
  par(mar = c(2.5,2.5,1,1))
  plot( nodes_mat[,2], nodes_mat[,1], pch = membership_true, col = membership_est, cex = 1, xlab='', ylab = '', xlim=c(0,6), ylim=c(0,1))
  
  
  # plot the circles
  # angel = seq(0, 2*pi, length.out=200)
  # x_center = centers[1,1]; y_center = centers[1,2]
  # points(y_center, x_center, pch = 21, bg='darkgray')
  # points( y_center+radius_thres1*sin(angel), x_center+radius_thres1*cos(angel), cex=0.1, col='darkgray')
}

# plot estimated mean pdf’s with true pdf's -------------------------------

n0_vec = rep(0, length(f_list))
for (i in 1:length(clusters)) {
  n0_vec[clusters[[i]]] = n0_mat[clusters[[i]], clusters[[i]][1]] - min(n0_mat[clusters[[i]], clusters[[i]][1]])
}
res = re_center_gd(f_list, clusters, n0_vec=n0_vec, 0.02, pp=pp)
f_center_list = res$f_center_list
n0_vec  = res$n0_vec

for (l in 1:length(clusters)) {
  f_center = f_center_list[[l]]

  if (pp)   {
    pdf_center = tail(f_center,length(t))
    plot(t, pdf_center, col = 'red', type='l', xlim = c(0,50), ylim=c(0,.4))
  }
  else   {
    pdf_center = obtain_pdf(tail(f_center,length(t)), t, bw=1)$density
    plot(pdf_center, col = 'red', xlim = c(0,50), ylim=c(0,.4), main = '')
  }
  
  for (i in (clusters[[l]])) {
    lines(t, shift(pdf_list[[i]], n0_vec[i], pad=0))
  }
  
  if (pp) lines(t, pdf_center, col = 'red', xlim = c(0,50), ylim=c(0,.4))
  else lines(pdf_center, col = 'red', xlim = c(0,50), ylim=c(0,.4))
}

  
#   -----------------------------------------------------------------------


# Plot node locations -----------------------------------------------------

case = 2

r = results2[[3]]
network = r$network
nodes_mat = network$nodes_mat


if (case==1)
{
  radius_thres = 1
  
  clus_size_1 = 4; clus_size_2 = 46
  centers = nodes_mat[1:clus_size_1,]
  
  dev.new(width=6,height=1.5,noRStudioGD = T)
  par(mar = c(2.5,2.5,1,1))
  plot( nodes_mat[,2], nodes_mat[,1], cex = .2, xlab='', ylab = '', xlim=c(0,6), ylim=c(0,1))
  points(nodes_mat[1:clus_size_1,2], nodes_mat[1:clus_size_1,1], col='red')
  
  # plot the circle
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points(y_center+radius_thres*sin(angel), x_center+radius_thres*cos(angel), cex=0.1, col='red')
}
if (case==2||case==3)
{
  radius_thres1 = 2
  radius_thres2 = 1
  
  clus_size_1 = 10; clus_size_2 = 10; clus_size_3 = 40
  centers = nodes_mat[1:(clus_size_1+clus_size_2),]
  
  dev.new(width=6,height=1.5,noRStudioGD = T)
  par(mar = c(2.5,2.5,1,1))
  plot( nodes_mat[,2], nodes_mat[,1], cex = .2, xlab='', ylab = '', xlim=c(0,6), ylim=c(0,1))
  points(nodes_mat[1:clus_size_1,2], nodes_mat[1:clus_size_1,1],  col='red')
  points(nodes_mat[1:clus_size_2+clus_size_1,2], nodes_mat[1:clus_size_2+clus_size_1,1], col='blue')
  
  # plot the circles
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points( y_center+radius_thres1*sin(angel), x_center+radius_thres1*cos(angel), cex=0.1, col='red')
  x_center = centers[1+clus_size_1,1]; y_center = centers[1+clus_size_1,2]
  points(y_center+radius_thres2*sin(angel), x_center+radius_thres2*cos(angel), cex=0.1, col='blue')
}



# plot the growing network ------------------------------------------------

case = 3
SEED = 2982
total_time = 50
t = seq(0, total_time, 0.05)

if (case==1) {
  network = generate_network(SEED, total_time)
}
if (case==2) {
  network = generate_network2(SEED, total_time)
}
if (case==3) {
  network = generate_network3(SEED, total_time)
}

nodes_mat = network$nodes_mat
edge_time_mat = network$edge_time_mat


{
  library(colorRamps)
  library(grDevices)
  library(fields)
  # colorbar = colorRamp(c(rgb(1,0,0,0), rgb(1,0,0,1)), alpha=T)
  colorbar = cm.colors(51)
  colorbar = blue2red(50)
  
  
  # dev.new(width=3, height=6, noRStudioGD = T)
  
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


# plot aligned cdf's & estimated mean cdf's by cluster -------------------------------

# total_time = 50
# t = seq(0, total_time, 0.05)
# 
case = 3
# r = results[[4]]

n0_vec = r$n0_ve
f_center_list = r$f_center_list
clusters = r$clusters
f_list = r$f_list

if (case==1)
{
    clus_col = c('red', 'black')
    colors = c(rep(clus_col[1], 4), rep(clus_col[2], 46))
    dev.new(width=6, height=4, noRStudioGD = T)
    par(mfrow=c(1,2))
    order = c(2,1)

    iter = 1
    for (clus_id in order) {
      col = clus_col[iter]
      plot(t, tail(shift(f_center_list[[clus_id]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col=col, xlab='time', ylab='f(t)')
      for (i in clusters[[clus_id]]) {
        lines(t, tail(shift(f_list[[i]], n0_vec[i]), length(t)), col=scales::alpha(colors[i],.3), lwd = .8)
      }
      lines(t, tail(shift(f_center_list[[clus_id]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col=col, xlab='time', ylab='f(t)')
      iter = iter+1
    }
}

if (case==2 || case==3)
{
    clus_col = c('red',  'blue', 'black')
    colors = c(rep(clus_col[1], 10), rep(clus_col[2], 10), rep(clus_col[3], 40))
    dev.new(width=8, height=4, noRStudioGD = T)
    par(mfrow=c(1,3))

    order = c(1,2,3)
    iter = 1
    for (clus_id in order) {
      col = clus_col[iter]
      plot(t, tail(shift(f_center_list[[clus_id]], 0, pp=pp), length(t)), type='l', lty='dashed', lwd=1.5,col=col, xlab='time', ylab='f(t)')
      for (i in clusters[[clus_id]]) {
        lines(t, tail(shift(f_list[[i]], n0_vec[i], pp=pp), length(t)), col=scales::alpha(colors[i],.3), lwd = .8)
      }
      lines(t, tail(shift(f_center_list[[clus_id]], 0, pp=pp), length(t)), type='l', lty='dashed', lwd=1.5,col=col, xlab='time', ylab='f(t)')
      iter = iter+1
    }
}



