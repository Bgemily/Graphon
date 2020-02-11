
rm(list=ls())
source('./cluster_point_proc.R')

# ARI ---------------------------------------------------------------------

get_ARI = function(membership_true, results, truncate)
{
  ARI = vector('numeric', length(results))
  for (i in 1:length(results)) {
    clusters_i = results[[i]]$clusters
    membership_i = unlist(clusters_i)
    for (j in 1:length(clusters_i)) {
      membership_i[clusters_i[[j]]] = j
    }
    ARI[i] = mclust::adjustedRandIndex(membership_i[1:truncate], membership_true[1:truncate])
  }
  return(ARI)
}

membership_true = c(rep(1,4), rep(2,46))
ARI1 = get_ARI(membership_true, results, 4)

membership_true2 = c(rep(1,4), rep(2,8), rep(3,46))
ARI2 = get_ARI(membership_true2, results2, 12)

membership_true3 = c(rep(1,4), rep(2,8), rep(3,46))
ARI3 = get_ARI(membership_true3, results3, 12)

ARI=ARI2
hist(ARI)
summary(ARI)
boxplot(ARI)

boxplot(c(ARI1, ARI2, ARI3)~c(rep(1,length(ARI1)), rep(2,length(ARI2)), rep(3,length(ARI3))))




# Initialization vs ARI ---------------------------------------------------

case = 1
if(case==1) {result = results; membership_true = c(rep(1,4), rep(2,46))}
if(case==2) {result = results2; membership_true = c(rep(1,4), rep(2,8), rep(3,46))}
if(case==3) {result = results3; membership_true = c(rep(1,4), rep(2,8), rep(3,46))}

ARI = get_ARI(membership_true, result, length(membership_true))
# ARI = get_ARI(membership_true, result, 12)

initials = t(sapply(c(1:length(result)), function(i)result[[i]][[6]]))
if (case==2||case==3) good_initials = apply(initials, 1, function(x)(min(x)<=4 && max(x)>12 && median(x)<=12 && median(x)>4))
if (case==1) good_initials = apply(initials, 1, function(x)(min(x)<=4 && max(x)>4))

par(mar = c(2.5, 2.5, 1, 1))
boxplot(ARI~factor(good_initials, levels=c("TRUE","FALSE"), labels=c("good",'fair')))

sum(good_initials==T); sum(good_initials==F); 

plotrix::multhist(list(ARI[good_initials], ARI[!good_initials]), freq=F)
legend('topleft',legend=c('good','fair'),col=c('gray20','gray90'),fill=c('gray30','gray90'))


# Plot node locations -----------------------------------------------------

case = 2

total_time = 50
t = seq(0, total_time, 0.05)
SEED = 203

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

if (case==1)
{
  radius_thres = 1
  
  clus_size_1 = 4; clus_size_2 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2))
  centers = cbind(c(rep(0.5,clus_size_1)), c(seq(0.8,5.2,length.out=clus_size_1)))
  
  dev.new(width=2,height=6,noRStudioGD = T)
  par(mar = c(2.5,2.5,1,1))
  plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '', xlim=c(0,1))
  points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  
  # plot the circle
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points(x_center+radius_thres*cos(angel), y_center+radius_thres*sin(angel), cex=0.1, col='red')
}
if (case==2||case==3)
{
  radius_thres1 = 2
  radius_thres2 = 1
  
  clus_size_1 = 4; clus_size_2 = 8; clus_size_3 = 46
  clusters_vec = c(rep(1,clus_size_1), rep(2,clus_size_2), rep(3, clus_size_3))
  centers = cbind(c(rep(0.5,clus_size_1), rep(c(0.7,0.3), clus_size_2/2)), c(seq(0.8,5.2,length.out=clus_size_1), seq(1,5,length.out=clus_size_2)))
  
  dev.new(width=2,height=6,noRStudioGD = T)
  par(mar = c(2.5,2.5,1,1))
  plot(nodes_mat[,1], nodes_mat[,2], cex = .2, xlab='', ylab = '', xlim=c(0,1))
  points(nodes_mat[1:clus_size_1,1], nodes_mat[1:clus_size_1,2], col='red')
  points(nodes_mat[1:clus_size_2+clus_size_1,1], nodes_mat[1:clus_size_2+clus_size_1,2], col='blue')
  
  # plot the circles
  angel = seq(0, 2*pi, length.out=200)
  x_center = centers[2,1]; y_center = centers[2,2]
  points(x_center+radius_thres1*cos(angel), y_center+radius_thres1*sin(angel), cex=0.1, col='red')
  x_center = centers[1+clus_size_1,1]; y_center = centers[1+clus_size_1,2]
  points(x_center+radius_thres2*cos(angel), y_center+radius_thres2*sin(angel), cex=0.1, col='blue')
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


# plot original cdf's -----------------------------------------------------

total_time = 50
t = seq(0, total_time, 0.05)

case = 1
r = results[[4]]
f_list = r$f_list
if (case==1)
{
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
}

if(case==2 || case==3)
{
    colors = c(rep('red', 4), rep('blue', 8), rep('black', 46))
    dev.new(width=8, height=4, noRStudioGD = T)
    par(mfrow=c(1,3))
    plot(t, tail(f_list[[1]], length(t)), type='l',col=colors[1], xlab='time', ylab='f(t)', main='Cluster 1')
    for (i in 2:4) {
      lines(t, tail(f_list[[i]], length(t)), col=colors[i])
    }
    plot(t, tail(f_list[[5]], length(t)), type='l',col=colors[5], xlab='time', ylab='f(t)', main='Cluster 2')
    for (i in 6:12) {
      lines(t, tail(f_list[[i]], length(t)), col=colors[i])
    }
    plot(t, tail(f_list[[13]], length(t)), type='l',col=colors[13], xlab='time', ylab='f(t)', main='Cluster 3')
    for (i in 14:58) {
      lines(t, tail(f_list[[i]], length(t)), col=colors[i])
    }
}


# plot aligned cdf's as well as the centers -------------------------------

total_time = 50
t = seq(0, total_time, 0.05)

case = 1
r = results[[4]]

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
    colors = c(rep(clus_col[1], 4), rep(clus_col[2], 8), rep(clus_col[3], 46))
    dev.new(width=8, height=4, noRStudioGD = T)
    par(mfrow=c(1,3))

    order = c(3,1,2)
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



# plot centers ------------------------------------------------------------

r = results3[[2]]
f_center_list = r$f_center_list

{
    # dev.new(width=6, height=6, noRStudioGD = T)
    # plot(t, tail(shift(f_center_list[[1]], 0), length(t)), ylim = c(0,1), type='l', lty='dashed', lwd=1.5, col='black', xlab='time', ylab='f(t)')
    # lines(t, tail(shift(f_center_list[[2]], 0), length(t)), type='l', lty='solid', lwd=1.5, col='red')
    # legend("topleft", legend=c("Cluster 1", "Cluster 2"), col=c("black", "red"), lty=c('dashed','dashed'), cex=0.8)
}

{
    dev.new(width=6, height=6, noRStudioGD = T)
    plot(t, tail(shift(f_center_list[[1]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='red', xlab='time', ylab='f(t)')
    lines(t, tail(shift(f_center_list[[2]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='black')
    lines(t, tail(shift(f_center_list[[3]], 0), length(t)), type='l', lty='dashed', lwd=1.5,col='blue')
    legend("topleft", legend=c("Cluster 1", "Cluster 2", 'Cluster 3'), col=c("red", "black", "blue"), lty=c('dashed','dashed','dashed'), cex=0.8)
}



##### plot of time delay #####

# This is not a good idea for the following reasons.
# First, the estimate of tau will be poor when misclassification appears.
# Second, we do not care about tau's (maybe).
# but we can shift min_tau to 0?

# Plot f_list -------------------------------------------------------------

case = 1
SEED = 203

total_time = 50
t = seq(0, total_time, 0.05)

if (case==1) network = generate_network(SEED, total_time)

if (case==2) network = generate_network2(SEED, total_time)

if (case==3) network = generate_network3(SEED, total_time)


edge_time_mat = network$edge_time_mat

f_list = get_emp_f_list(edge_time_mat, t)
for (i in c(1:10)) {
  plot(f_list[[i]], type='l', main=i)
}





