library(data.table)
library(R.matlab)


path.list=list.files('../Zebrafish_spinal_cord_development/FunctionalData/');
k=1;
path=path.list[[k]]
dat=readMat(paste('../Zebrafish_spinal_cord_development/FunctionalData/',path,'/profile.mat',sep=''))

locs.all=cbind(dat$x,dat$y,dat$z)
rm(dat)

edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
edge.time=edge.time[,-1]


avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
avai.inds=avai.inds[,-1];

member.ship = as.matrix(read.csv(paste('./FunctionalData/',path,'/MembShip.csv',sep='')))
member.ship=member.ship[,-1];


locs=locs.all[avai.inds,]

color.list=c('red','blue','green') # Add more colors if you want to

#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
scatterplot3d(x=locs[,1], y=locs[,2], z=locs[,3],xlim=range(locs[,1]),ylim=range(locs[,2]),zlim=range(locs[,3]),color=color.list[member.ship],pch=16)

for (t in seq(120,80,-5)) {
  plot(x=locs[,1], y=locs[,2], col=member.ship+1,pch=16, xlab = t)
  for (i in 1:dim(locs)[1]) {
    for (j in 1:dim(locs)[1]) {
      if (edge.time[i, j]<t+5 & edge.time[i,j]>t) {
        # if (member.ship[i]==2 & member.ship[j]==2 & edge.time[i,j]< Inf) {
        lines(x = c(locs[i,1], locs[j,1]), y = c(locs[i,2], locs[j,2]))
      }
    }
  }
}

idx = c(13,18,21,23,24,26,29)
points(x=locs[idx,1], y=locs[idx,2], col=1,pch=12, xlab = t)

first_edge_time = matrix(Inf, nrow=nrow(edge.time), ncol=ncol(edge.time))
for (i in 1:nrow(edge.time)) {
  idx_tmp = order(edge.time[i,])[2] # the second smallest element 
  first_edge_time[i, idx_tmp] = edge.time[i,idx_tmp]
  first_edge_time[idx_tmp,i] = edge.time[i,idx_tmp]
}

for (t in seq(280,80,-10)) {
  plot(x=locs[,1], y=locs[,2], col=member.ship+1,pch=16, xlab = t)
  for (i in 1:dim(locs)[1]) {
    for (j in 1:dim(locs)[1]) {
      if (first_edge_time[i, j]<t+10 & first_edge_time[i,j]>=2) {
        # if (member.ship[i]==2 & member.ship[j]==2 & edge.time[i,j]< Inf) {
        lines(x = c(locs[i,1], locs[j,1]), y = c(locs[i,2], locs[j,2]))
      }
    }
  }
}


plot(0,xlim=c(0,340),ylim=c(0,0.005))
for(i in 2:dim(edge.time)[1]){
  if(length(which(first_edge_time[i,-i]<Inf))>=2)
    lines(density(first_edge_time[i,-i]))
}


