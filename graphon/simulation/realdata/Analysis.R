#### Basic connectivity inference
library(data.table)
path.list=list.files('../processed_FunctionalData/');


conn.interval=200; # =50 seconds at 4 Hz
rho=0.5; # an arbitrary threshold...

window_length = 240 # = 1min at 4 Hz
window_step = 240 # = 1min at 4Hz

# Calculate the covariance matrix 
# Take the inverse of the covariance matrix 

path.list=list.files('../processed_FunctionalData/');

for(k in 1:length(path.list)){
  path=path.list[[k]]
  dat.dFF=as.matrix(fread(paste('../processed_FunctionalData/',path,'/dFF.csv',sep='')))
  dat.dFF=dat.dFF[,-1]
  # Simply calculate the covariance among dFF within the time periods
  # Obtain the list of missing traces 
  na.inds= is.na(dat.dFF[,1]);
  reduced.dFF=dat.dFF[!na.inds,];
  
  # Calculate the covariance/correlation matrices for each time interval
  
  # n.intervals= round(dim(reduced.dFF)[2]/conn.interval,digits=0)
  # interval.list=matrix(1:dim(reduced.dFF)[2],ncol=conn.interval,byrow=T)
  # interval.list=interval.list[1:n.intervals,];
  
  n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
  interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
  interval.list = interval.list + (1:n.intervals-1)*window_step
  
  
  cor.full=array(dim=c(n.intervals,dim(reduced.dFF)[1],dim(reduced.dFF)[1]))
  for(i in 1:n.intervals){
    dat.temp=reduced.dFF[,interval.list[i,]];
    cor.full[i,,]=cor(x=t(dat.temp));  
  }
  
  # Apply glasso on each one of the covariance matrix 
  # Too slow
  # library(glasso)
  # cov.tmp=cov.full[200,,]
  # glassopath(s=cov.tmp)
  
  # Thresholding the correlation matrix to obtain the connectivity matrix 
  # nodes are considered as connected if their mean correlation in next 5 min is larger than the threshold
  adj.full=cor.full;
  for(i in 1:(n.intervals-4)){
    adj.full[i,,]=apply(cor.full[i:(min(i+4,n.intervals)),,,drop=F], c(2,3), mean)>rho;
  }
  
  # Find the minimum connecting time for the neurons 
  edge.time=matrix(Inf,nrow=dim(reduced.dFF)[1],ncol=dim(reduced.dFF)[1])
  for(i in 1:dim(reduced.dFF)[1]){
    for(j in 1:dim(reduced.dFF)[1]){
      edge.time[i,j]=min(which(adj.full[,i,j]==1))
    }
  }
  
  
  write.csv(which(!na.inds),paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep=''),col.names=F)
  
  write.csv(edge.time,paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep=''),col.names=F)
  
}

## Load the edge time matrix
path.list=list.files('./FunctionalData/');
path=path.list[[1]];
edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
edge.time=edge.time[,-1]

image(edge.time)

bw=10
plot(density(edge.time[1,-1]),xlim=c(0,340),ylim=c(0,0.02))
for(i in 2:dim(edge.time)[1]){
  if(length(which(edge.time[i,-i]<Inf))>=2)
    lines(density(edge.time[i,-i]))
}
i=16;lines(density(edge.time[i,-i],bw=5), col=2,xlim=c(0,340),ylim=c(0,0.02))