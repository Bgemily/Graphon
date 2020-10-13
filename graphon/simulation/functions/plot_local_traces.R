
##### plot local traces

library(ggridges)
library(ggplot2)
library(reshape2)
library(gridExtra)

plot_local_traces = function(locs, reduced.dFF, edge_time_mat=NULL, snapshot_mins=c(30,120,240), 
                             window_length=1, id_vec=NULL, show_plot=FALSE,  scale=2){
  order.temp = order(locs[,1], decreasing = TRUE)
  order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])
  
  if(!is.null(edge_time_mat)){
    tmp_ind = 1:nrow(locs)
    min_edge_time = sapply(1:length(tmp_ind), function(i)min(edge_time_mat[tmp_ind,tmp_ind][i,-i]))
    order.temp = order(min_edge_time); 
    order.temp = rev(order.temp)
    order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])
  }
  
  
  
  if(!is.null(id_vec))
    id_vec = sapply(id_vec, function(id)which(order.temp==id))
  
  timefield = lapply(snapshot_mins, function(min)min*240+1:(240*window_length))
  g.temp = list()
  for (t in timefield) {
    mat.temp = reduced.dFF[,t]
    mat.temp = mat.temp[order.temp,]
    LR = (locs[order.temp,2]<0)
    thick = vector(length=nrow(reduced.dFF))
    thick[id_vec] = TRUE
    color = rep("0",nrow(reduced.dFF))
    if (!is.null(id_vec)){
      color[id_vec[1]] = "1"; color[id_vec[2]] = "2"; color[id_vec[3]] = "3"
    }
    df.temp = data.frame(id=1:nrow(mat.temp), mat.temp)
    df.temp = reshape2::melt(df.temp, id = "id")
    g<-ggplot(df.temp, aes(x = variable, y = id, height = value, group=id, fill=LR[id], size=thick[id], color=color[id])) +
      geom_ridgeline(scale=scale, min_height=-0.1) +
      scale_size_manual(values = c("TRUE"=1, "FALSE"=0.1))+
      scale_color_manual(values = c("0"="gray","1"=viridis::viridis(3)[1],"2"=viridis::viridis(3)[2],"3"=viridis::viridis(3)[3]))+
      xlab(paste(round((t[1]-1)/240/60,digits = 1), "h"))+
      theme(legend.position = "none", axis.text.x=element_blank(), 
            axis.ticks.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank(), axis.title.y = element_blank()) 
    g.temp = c(g.temp, list(g))
  }
  if(show_plot)
    do.call("grid.arrange", c(g.temp, ncol=length(g.temp)))
  else
    return(g.temp)
}
