
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
clusters_history = clus_result$clusters_history
res_overclus = clus_result$res_overclus
n0_vec = clus_result$n0_vec



membership_est = clus2mem(clusters_est)


par(mar=c(2.5,2.5,.5,.5))

# ARI ---------------------------------------------------------------------


results = results2
membership_true = results[[1]]$network$membership_true

clusters_list = lapply(results, function(x)x$clus_result$clusters)
ARI = get_ARI(membership_true, clusters_list)

clusters_list_exaclus = lapply(results, function(x)x$clus_result$clusters_history[[1]])
ARI_exaclus = get_ARI(membership_true, clusters_list_exaclus)

clusters_list_iter_1 = lapply(results, function(x)x$clus_result$clusters_history[[3]])
ARI_iter_1 = get_ARI(membership_true, clusters_list_iter_1)
clusters_list_iter_2 = lapply(results, function(x)x$clus_result$clusters_history[[5]])
ARI_iter_2 = get_ARI(membership_true, clusters_list_iter_2)
clusters_list_iter_3 = lapply(results, function(x)x$clus_result$clusters_history[[7]])
ARI_iter_3 = get_ARI(membership_true, clusters_list_iter_3)





plot_jitter_boxplot = function(data)
{  
  library(ggplot2)
  library(dplyr)
  library(viridis)
  data = reshape2::melt(data, id.vars=NULL, variable.name="name")
  sample_size = data %>% group_by(name) %>% summarize(num=n())
  data %>%
    left_join(sample_size) %>%
    # mutate(myaxis = paste0(name)) %>%
    ggplot( aes(x=name, y=value, fill=name)) +
    geom_violin(width=1, alpha=0.3) +
    geom_boxplot(width=0.3, alpha=0.7) +
    scale_fill_viridis(discrete = TRUE) +
    ylim(c(0,1))+
    theme_light() +
    theme(
      axis.title = element_blank(),
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("") 
    # coord_flip()
  
}


plot_jitter_boxplot(data=data.frame(ARI_iter_1, ARI_iter_2, ARI_iter_3 , ARI))
plot_jitter_boxplot(data=data.frame(ARI_iter_2))

plot_jitter_boxplot(data=data.frame(ARI))


# plot_jitter_boxplot(data=data.frame(ARI_exaclus, ARI_overclus_1, ARI_overclus_2, ARI_overclus_3 , ARI))
# 
# 
# plot_jitter_boxplot(data=data.frame(ARI_overclus_3.3, ARI_overclus_3.4, ARI_overclus_3.5, ARI_overclus_3.6))



# Plot estimated connecting patterns -------------------------------------

results = results2


pdf_true_array = fun2pdfarray(true_pdf_fun_list = results[[1]]$network$true_pdf_fun_list,
                               tau_mat = results[[1]]$network$tau_mat, 
                              membership_true = results[[1]]$network$membership_true,
                              t_vec = r$network$t_vec)
# times pdf_true_array by connecting probabilities
r = results[[2]]
edge_time_mat = r$network$edge_time_mat
clus_degree = get_clus_degree_mat(edge_time_mat = edge_time_mat, clusters = mem2clus(r$network$membership_true))
for (q in 1:dim(pdf_true_array)[1]) {
  for (l in 1:dim(pdf_true_array)[2]) {
    pdf_true_array[q,l, ] = pdf_true_array[q,l, ] * clus_degree[q,l] / (sum(r$network$membership_true==q)*sum(r$network$membership_true==l)*(1-0.5*I(q==l)))
  }
}


pdf_array_list = lapply(results, function(r)r$clus_result$center_pdf_array)
clusters_list = lapply(results, function(x)x$clus_result$clusters)

# permutate clusters and pdf_array's for each subject
index = 1:10
res = match_clusters(clusters_list = clusters_list[index], pdf_array_list = pdf_array_list[index],
                     pdf_true_array = pdf_true_array)
clusters_list = res$clusters_list
pdf_array_list = res$pdf_array_list




r = results[[1]]
plot_pdf_array(pdf_array_list = pdf_array_list[index], pdf_true_array = pdf_true_array, 
               t_vec = r$network$t_vec, y_lim = c(0,0.25)) 
plot_pdf_array(pdf_array_list = center_pdf_array, pdf_true_array = pdf_true_array, 
               t_vec = r$network$t_vec, y_lim = c(0,0.25)) 





# Plot estimated connecting patterns (real data) -------------------------------------

results = list(result)


pdf_array_list = lapply(results, function(r)r$clus_result$center_pdf_array)


r = results[[1]]
plot_pdf_array(pdf_array_list = pdf_array_list, pdf_true_array = NULL, 
               t_vec = r$network$t_vec, y_lim = c(0,0.02)) 


write.csv(clus2mem(result$clus_result$clusters), paste('../processed_FunctionalData/',path,'/MembShip.csv',sep=''),col.names=F)



# Over-clustering -------------------------------------------------

r = results2[[29]]
edge_time_mat = r$network$edge_time_mat
t_vec = r$network$t_vec
bw = 1
membership_true = r$network$membership_true

clusters_exaclus = r$clus_result$clusters_history[[1]]
membership_exaclus = clus2mem(clusters_exaclus)

clusters_iter_1 = r$clus_result$clusters_history[[2]]
membership_iter_1 = clus2mem(clusters_iter_1)
clusters_merged_1 = r$clus_result$clusters_history[[3]]
membership_merged_1 = clus2mem(clusters_merged_1)
clusters_iter_2 = r$clus_result$clusters_history[[4]]
membership_iter_2 = clus2mem(clusters_iter_2)
clusters_merged_2 = r$clus_result$clusters_history[[5]]
membership_merged_2 = clus2mem(clusters_merged_2)
clusters_iter_3 = r$clus_result$clusters_history[[6]]
membership_iter_3 = clus2mem(clusters_iter_3)
clusters_iter_4 = r$clus_result$clusters_history[[8]]
membership_iter_4 = clus2mem(clusters_iter_4)


plot_tsne = function(dist_mat, membership_true, membership_est){
  set.seed(81)
  tsne = Rtsne::Rtsne(dist_mat, is_distance=TRUE, perplexity=20)$Y
  # plot(tsne, pch=membership_true, col=membership_iter_1)
  df = data.frame(x=tsne[,1], y=tsne[,2], 
                  true_clus = as.factor(membership_true), 
                  est_clus = as.factor(membership_est))
  library(ggplot2)
  ggplot(df, aes(x=x, y=y, color=est_clus, shape=true_clus)) + 
    geom_point(size=2) +
    scale_shape_manual(values=c(1,2,4)) + 
    theme( axis.title = element_blank() ,
           axis.ticks = element_blank(),
           axis.text = element_blank())
  
}



N_node = nrow(edge_time_mat)
node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                    n0_vec = numeric(N_node), t_vec = t_vec, bw = bw)
n0_vec_init = est_n0_vec(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), t_vec = t_vec, bw = bw)
aligned_pdf_mat = t(sapply(1:N_node, function(i)shift(node_pdf_array[i,1,], n0_vec_init[i], pp=TRUE)))
dist_mat = rdist::pdist(aligned_pdf_mat)
plot_tsne(dist_mat = dist_mat,  membership_true = membership_true, 
          membership_est = membership_iter_1)
get_one_ARI(memb_est_vec = membership_iter_1,  memb_true_vec = membership_true)



get_pairwise_dist = function(edge_time_mat, clusters, t_vec=seq(0, 50, 0.05), bw=1)
{
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                       n0_vec = n0_vec, t_vec = t_vec, bw = bw)
  degree_mat = get_node_degree_mat(edge_time_mat = edge_time_mat, clusters = clusters)
  dist_mat = pairwise_dist_mat(pdf_array = node_pdf_array, degree_mat = degree_mat, t_unit = 0.05)$dist_mat
  return(dist_mat)
}

dist_mat_iter_1 =  get_pairwise_dist(edge_time_mat = edge_time_mat, clusters = clusters_iter_1)
plot_tsne(dist_mat_iter_1, membership_true, membership_iter_2)
plot_tsne(dist_mat_iter_1, membership_true, membership_merged_2)



dist_mat_iter_2 =  get_pairwise_dist(edge_time_mat = edge_time_mat, clusters = clusters_iter_2)
plot_tsne(dist_mat_iter_2, membership_true, membership_iter_3)


dist_mat_iter_3 =  get_pairwise_dist(edge_time_mat = edge_time_mat, clusters = clusters_iter_3)
plot_tsne(dist_mat_iter_3, membership_true, membership_iter_4)



plot_conn_patt = function(edge_time_mat, clusters, t_vec=seq(0, 50, 0.05), bw=1)
{
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  center_pdf_array = get_center_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                      n0_vec = n0_vec, t_vec = t_vec, bw = bw)
  plot_pdf_array(center_pdf_array)
}
plot_conn_patt(edge_time_mat = edge_time_mat, clusters = clusters_iter_2)



# case1: pdf vs cdf --------------------------------------------------------------


tmp.t = seq(0, 30,0.05)
tmp.tau = tmp.t - 10
tmp.f1 = sapply(tmp.t, dnorm, mean=10, sd=1)
tmp.f2 = sapply(tmp.t, dnorm, mean=20, sd=1)

theta1 = fft(tmp.f1); theta2 = fft(tmp.f2)

tmp.dist = sapply(tmp.tau/0.05, function(n0)distance(theta2, theta1, n0, pp=TRUE))

tmp.df = data.frame(t=tmp.t, tau=tmp.tau, f1=tmp.f1, f2=tmp.f2, dist=tmp.dist^2)

library(ggplot2)
p1 = ggplot(tmp.df) + 
      geom_line(aes(x=t, y=f1), linetype=2) + 
      geom_line(aes(x=t, y=f2) ) + 
      xlab("x")+ylab("f(x)")
p2 = ggplot(tmp.df) + 
      geom_line(aes(x=tau, y=dist))+
      xlab("shift")+ylab(expression(d^2))

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))



tmp.cdf1 = cumsum(tmp.f1)/sum(tmp.f1)
tmp.cdf2 = cumsum(tmp.f2)/sum(tmp.f2)

cdf.theta1 = fft(tmp.cdf1)
cdf.theta2 = fft(tmp.cdf2)
N = length(cdf.theta1)
k = seq(1, N-1)
  {
    cdf.theta1_prime = c(cdf.theta1[1], cdf.theta1[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
    cdf.theta2_prime = c(cdf.theta2[1], cdf.theta2[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  }
  

cdf.tmp.dist = sapply(tmp.tau/0.05, function(n0)distance(cdf.theta2_prime, cdf.theta1_prime, n0, pp=FALSE))

cdf.tmp.df = data.frame(t=tmp.t, tau=tmp.tau, f1=tmp.cdf1, f2=tmp.cdf2, dist=cdf.tmp.dist^2)

library(ggplot2)
p1 = ggplot(cdf.tmp.df) + 
  geom_line(aes(x=t, y=f1), linetype=2) + 
  geom_line(aes(x=t, y=f2)) + 
  xlab("x")+ylab("F(x)")
p2 = ggplot(cdf.tmp.df) + 
  geom_line(aes(x=tau, y=dist), )+
  xlab("shift")+ylab(expression(d^2))

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))



# case2: pdf    ###########################################################


tmp.t = seq(0, 50,0.05)
tmp.tau = tmp.t - 10
tmp.f1 = sapply(tmp.t, function(x) dnorm(x, 10, 1))
tmp.f2 = sapply(tmp.t, function(x) dunif(x, 5,15))
tmp.f3 = sapply(tmp.t, function(x) dnorm(x, 10, 1)*0.8 + dnorm(x,40,2)*0.2)


theta1 = fft(c(rep(0,length(tmp.f1)),tmp.f1)); 
theta2 = fft(c(rep(0,length(tmp.f2)),tmp.f2))
theta3 = fft(c(rep(0,length(tmp.f3)),tmp.f3))

tmp.dist12 = sapply(tmp.tau/0.05, function(n0)distance(theta2, theta1, n0, pp=TRUE))
tmp.dist13 = sapply(tmp.tau/0.05, function(n0)distance(theta3, theta1, n0, pp=TRUE))

tmp.df = data.frame(t=tmp.t, tau=tmp.tau, f1=tmp.f1, f2=tmp.f2, f3=tmp.f3, dist12=tmp.dist12^2, dist13=tmp.dist13^2)

library(ggplot2)
p1 = ggplot(tmp.df) + 
  geom_line(aes(x=t, y=f1), linetype=2) + 
  geom_line(aes(x=t, y=f2), color=2 ) + 
  geom_line(aes(x=t, y=f3), color=3 ) + 
  xlab("x")+ylab("f(x)")
p2 = ggplot(tmp.df) + 
  geom_line(aes(x=tau, y=dist12), color=2)+
  geom_line(aes(x=tau, y=dist13), color=3)+
  xlab("shift")+ylab(expression(d^2))

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))


# case2: cdf -----------------------------------------------------------

tmp.cdf1 = cumsum(tmp.f1)/sum(tmp.f1)
tmp.cdf2 = cumsum(tmp.f2)/sum(tmp.f2)
tmp.cdf3 = cumsum(tmp.f3)/sum(tmp.f3)

cdf.theta1 = fft(c(rep(0,length(tmp.cdf1)),tmp.cdf1))
cdf.theta2 = fft(c(rep(0,length(tmp.cdf2)),tmp.cdf2))
cdf.theta3 = fft(c(rep(0,length(tmp.cdf3)),tmp.cdf3))
N = length(cdf.theta1)
k = seq(1, N-1)
{
  cdf.theta1_prime = c(cdf.theta1[1], cdf.theta1[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  cdf.theta2_prime = c(cdf.theta2[1], cdf.theta2[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
  cdf.theta3_prime = c(cdf.theta3[1], cdf.theta3[2:N]+exp(-1i*2*pi*k)/(1-exp(-1i*2*pi*k/N)))
}


cdf.tmp.dist12 = sapply(tmp.tau/0.05, function(n0)distance(cdf.theta2_prime, cdf.theta1_prime, n0, pp=FALSE))
cdf.tmp.dist13 = sapply(tmp.tau/0.05, function(n0)distance(cdf.theta3_prime, cdf.theta1_prime, n0, pp=FALSE)) 

cdf.tmp.df = data.frame(t=tmp.t, tau=tmp.tau, f1=tmp.cdf1, f2=tmp.cdf2, f3=tmp.cdf3, dist12=cdf.tmp.dist12^2, dist13=cdf.tmp.dist13^2)

library(ggplot2)
p1 = ggplot(cdf.tmp.df) + 
  geom_line(aes(x=t, y=f1), linetype=2) + 
  geom_line(aes(x=t, y=f2), col=2) + 
  geom_line(aes(x=t, y=f3), col=3) + 
  xlab("x")+ylab("F(x)")
p2 = ggplot(cdf.tmp.df) + 
  geom_line(aes(x=tau, y=dist12), col=2)+
  geom_line(aes(x=tau, y=dist13), col=3)+
  xlab("shift")+ylab(expression(d^2)) 

p2.zoom = p2 + coord_cartesian(xlim = c(-2,2), ylim = c(0,2)) +
          theme(legend.position="none", 
                # axis.line=element_blank(),
                # axis.text.x=element_blank(), 
                axis.text.y=element_blank(),
                # axis.ticks=element_blank(),
                axis.title.x=element_blank(),axis.title.y=element_blank(),
                # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_rect(color='red', fill="white"),
                # plot.margin = unit(c(0,0,-6,-6),"mm")
                )
          
p3 = p2 + annotation_custom(ggplotGrob(p2.zoom), xmin=-10, xmax=10, ymin=15, ymax=40)


library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p3), size = "last"))





# Plot nodes locations with cluster results -------------------------------

library(R.matlab)

path.list=list.files('/Users/bgemily/Documents/Academic/SC/graphon/FunctionalData/');

for(k in 1:length(path.list)){
  path=path.list[[k]]
  dat<- readMat(paste('/Users/bgemily/Documents/Academic/SC/graphon/FunctionalData/',
                      path,'/profile.mat',sep=''));
  
  dat.activity=dat$profile.all;
  n.neurons = dim(dat.activity)[1]; 
  n.timeframes=dim(dat.activity)[2];
  
  write.csv(dat.dFF,paste('./FunctionalData/',path,'/dFF.csv',sep=''),col.names=F)
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





