
rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Setup  ------------------------------------------------------------------

total_time = 50
time_unit = 0.05
t = seq(0, total_time, 0.05)
pp = TRUE

r = results3[[1]]

f_list = r$f_list
f_center_list = r$f_center_list
clusters = r$clusters
# clusters = r$clusters_old
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


plot_jitter_boxplot(ARI=c(ARI3_old, ARI3),group=c(rep('',length(ARI3_old)),rep('cont.cluster',length(ARI3))))
plot_jitter_boxplot(ARI3_cont)

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


# plot center_pdf_array ------------------------------------------


plot(cmdscale(dist_mat), pch=membership_true, col=membership_est)
plot_pattern_matrix(clusters, edge_time_mat, n0_mat) # plot_center_pdf_array


# over-cluster
k=6
W = exp(-dist_mat^2/median(dist_mat)^2)
membership_overclus = spectral_clustering(W, k)
plot(cmdscale(dist_mat), pch=membership_true, col=membership_overclus)

clusters_overclus = list()
for (i in 1:k) {
  clusters_overclus[[i]] = which(membership_overclus==i)
}

# merge clusters, and plot clustering result in 2D plot.
clusters_merged = merge_clusters(clusters_overclus, edge_time_mat, n0_mat, 3)

membership_merged = unlist(clusters_merged)
for (j in 1:length(clusters_merged)) {
  membership_merged[clusters_merged[[j]]] = j
}
plot(cmdscale(dist_mat), pch=membership_true, col=membership_merged)


plot_pattern_matrix(clusters_overclus, edge_time_mat, n0_mat)
plot_pattern_matrix(clusters_merged, edge_time_mat, n0_mat)


# compare exact-clustering results and over-clustering results ---------------------------------------------------

k_trueclus = 3
k_overclus = 6
results = results3
ARI = ARI3

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

