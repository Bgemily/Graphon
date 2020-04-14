
rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Setup  ------------------------------------------------------------------

total_time = 50
time_unit = 0.05
t = seq(0, total_time, 0.05)

r = results2[[1]]

network = r$network
edge_time_mat = network$edge_time_mat
node_loc_mat = network$node_loc_mat
tau_vec = network$tau_vec
tau_mat = network$tau_mat
true_pdf_fun_list = network$true_pdf_fun_list
membership_true = network$membership_true
t_vec = network$t_vec
dist_thres = network$dist_thres
pairwise_dist = network$pairwise_dist

clusters_true = mem2clus(membership_true)


clus_result = r$clus_result
clusters_est = clus_result$clusters
center_pdf_array = clus_result$center_pdf_array
clusters_specc = clus_result$clusters_specc
clusters_specc_exaclus = clus_result$clusters_specc_exaclus

membership_est = clus2mem(clusters_est)
membership_specc = clus2mem(clusters_specc)
membership_specc_exaclus = clus2mem(clusters_specc_exaclus)



par(mar=c(2.5,2.5,.5,.5))

# ARI ---------------------------------------------------------------------

results = results2

membership_true = results[[1]]$network$membership_true
clusters_list = lapply(results, function(x)x$clus_result$clusters)
ARI = get_ARI(membership_true, clusters_list)

clusters_list_specc = lapply(results, function(x)x$clus_result$clusters_specc)
ARI_specc = get_ARI(membership_true, clusters_list_specc)

clusters_list_specc_exaclus = lapply(results, function(x)x$clus_result$clusters_specc_exaclus)
ARI_specc_exact = get_ARI(membership_true, clusters_list_specc_exaclus)

plot_jitter_boxplot(ARI=c(ARI_specc, ARI_specc_exact, ARI),
                    group=c(rep('spec (k=5) and merge',length(ARI_specc)),
                            rep('spec (k=3)',length(ARI_specc_exact)), 
                            rep('spec and k-means',length(ARI))) )




# Show estimated center_pdf_array (for one subject) ------------------------------------------

# need to find permutation of clusters that matches the est_pdf_array and true_pdf_array

pdf_true_array = fun2pdfarray(true_pdf_fun_list, tau_mat, membership_true)
plot_pdf_array(center_pdf_array, pdf_true_array)



# Confidence band of center_pdf_array -------------------------------------

results = results2
pdf_true_array = fun2pdfarray(results[[1]]$network$true_pdf_fun_list, 
                              results[[1]]$network$tau_mat, results[[1]]$network$membership_true)
pdf_array_list = lapply(results, function(r)r$clus_result$center_pdf_array)
clusters_list = lapply(results, function(x)x$clus_result$clusters)

# permutate clusters and pdf_array's for each subject
res = match_clusters(clusters_list = clusters_list, pdf_array_list = pdf_array_list, 
                     pdf_true_array = pdf_true_array)
clusters_list = res$clusters_list
pdf_array_list = res$pdf_array_list


alpha = 0.05
N_clus = dim(pdf_true_array)[1]
pdf_upper_array = array(dim = dim(pdf_true_array))
pdf_lower_array = array(dim = dim(pdf_true_array))
for (q in 1:N_clus) {
  for (l in 1:N_clus) {
    ql_pdf_mat = sapply(pdf_array_list, function(x)x[q,l,]) # len(t_vec)*n
    pdf_upper_array[q,l,] = apply(ql_pdf_mat, 1, quantile, 1-alpha/2)
    pdf_lower_array[q,l,] = apply(ql_pdf_mat, 1, quantile, alpha/2)
  }
}

plot_pdf_array(pdf_upper_array, pdf_true_array = pdf_true_array, pdf_array_2 = pdf_lower_array)



# Explain over-clustering -------------------------------------------------

par(mar=c(0,0,0,0))

X = cmdscale(dist_mat, k=4)
tmp = tsne::tsne(X, k=2)
tsne.2->tmp

plot(tmp)
plot(tmp, pch=membership_true, col=membership_specc_exaclus)
plot(tmp, pch=membership_true, col=membership_overclus)
plot(tmp, pch=membership_true, col=membership_specc)

tmp->tsne.2




# Explain node types --------------------------------------------------





# Plot nodes locations with cluster results -------------------------------



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


