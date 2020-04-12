plot_jitter_boxplot = function(ARI, group=1)
{  
  library(ggplot2)
  data_frame = data.frame(ARI, group=group)
  ggplot(data_frame, aes(x=group, y=ARI, group=group)) + 
    geom_boxplot(outlier.shape=NA)+
    # geom_violin()
    geom_jitter(position=position_jitter(width=.1, height=0))
}
