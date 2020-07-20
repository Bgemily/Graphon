library(data.table)
library(R.matlab)


path.list=list.files('../Zebrafish_spinal_cord_development/FunctionalData/');
k=8;
path=path.list[[k]]
dat=readMat(paste('../Zebrafish_spinal_cord_development/FunctionalData/',path,'/profile.mat',sep=''))

locs.all=cbind(dat$x,dat$y,dat$z)
mnx = dat$mnx
rm(dat)

edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
edge.time=edge.time[,-1]


avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
avai.inds=avai.inds[,-1];

# member.ship = as.matrix(read.csv(paste('./FunctionalData/',path,'/MembShip.csv',sep='')))
# member.ship=member.ship[,-1];

dat.dFF=as.matrix(fread(paste('../processed_FunctionalData/',path,'/dFF.csv',sep='')))
dat.dFF=dat.dFF[,-1]
reduced.dFF=dat.dFF[avai.inds,];

locs=locs.all[avai.inds,]

# color.list=c('red','blue','green') # Add more colors if you want to
# 
# #install.packages("scatterplot3d") # Install
# library("scatterplot3d") # load
# scatterplot3d(x=locs[,1], y=locs[,2], z=locs[,3],xlim=range(locs[,1]),ylim=range(locs[,2]),zlim=range(locs[,3]),color=color.list[member.ship],pch=16)
# 
# 
# for (t in seq(120,80,-5)) {
#   plot(x=locs[,1], y=locs[,2], col=mnx+1,pch=16, xlab = t)
#   for (i in 1:dim(locs)[1]) {
#     for (j in 1:dim(locs)[1]) {
#       if (edge.time[i, j]<t+5 & edge.time[i,j]>1) {
#         # if (member.ship[i]==2 & member.ship[j]==2 & edge.time[i,j]< Inf) {
#         lines(x = c(locs[i,1], locs[j,1]), y = c(locs[i,2], locs[j,2]))
#       }
#     }
#   }
# }
# 
# plot(0,xlim=c(0,340),ylim=c(0,0.05))
# for(i in which(mnx==0)){
#   if(length(which(edge.time[i,-i]<Inf))>=2)
#     lines(density(edge.time[i,-i]), col=mnx[i]+1)
# }
# 


##### ridgelineplot

library(ggridges)
library(ggplot2)
library(reshape2)

# 20min: 4800 
order.temp = order(locs[,1])
order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])

timefield = list(1600+1:240, 19000+1:240, 58000+1:240)
g.temp = list()
for (t in timefield) {
  mat.temp = reduced.dFF[,t]
  mat.temp = mat.temp[order.temp,]
  LR=(locs[order.temp,2]<0)
  df.temp = data.frame(id=1:nrow(mat.temp), mat.temp)
  df.temp = melt(df.temp, id = "id")
  g<-ggplot(df.temp, aes(x = variable, y = id, height = value, group=id, fill=LR[id])) +
    geom_ridgeline(scale=5, size=0.3, min_height=-0.1) +
    theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks.x = element_blank()) 
  g.temp = c(g.temp, list(g))
}
library(gridExtra)
grid.arrange(g.temp[[1]], g.temp[[2]], g.temp[[3]], ncol = 3)


## plot network
library(igraph)
adj.tmp = edge.time<=(68000/240+1) & edge.time>1
adj.tmp = adj.tmp[order.temp, order.temp]
colnames(adj.tmp) = rownames(adj.tmp) = 1:nrow(edge.time)
network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected")
plot(network.tmp, vertex.label=NA, palette=scales::hue_pal()(2), asp=.5, layout=locs[order.temp,], vertex.size=3,  vertex.color=(locs[order.temp,2]<0)+1)



## Load the edge time matrix

image(edge.time)

bw=10
edge.time.tmp = edge.time
i=2; plot(density(edge.time.tmp[i,-i]),xlim=c(0,340),ylim=c(0,0.04))
for(i in 2:nrow(edge.time.tmp)){
  if(length(which(edge.time.tmp[i,-i]<Inf))>=2)
    lines(density(edge.time.tmp[i,-i]))
}
i=order.temp[20];lines(density(edge.time.tmp[i,-i],bw=5), col=2,xlim=c(0,340),ylim=c(0,0.02))


act.time = edge.time[,1]
for (i in 1:nrow(edge.time)) {
  act.time[i] = min(edge.time[i,-i])
}
hist(act.time[mnx==0])
