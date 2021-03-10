
# Visualize summary statistics (median & quantile)  ---------------------------------------------------------------------


plot_jitter_boxplot = function(data, ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)), 
                               ylab=NULL, xlab=NULL){  
  library(ggplot2)
  library(dplyr)
  library(viridis)
  ylim = ylim
  data = reshape2::melt(data, id.vars=NULL, variable.name="name")
  sample_size = data %>% group_by(name) %>% summarize(num=n())
  data %>%
    left_join(sample_size) %>%
    # mutate(myaxis = paste0(name)) %>%
    ggplot( aes(x=name, y=value, fill=name)) +
    geom_violin(width=1, alpha=0.3) +
    geom_boxplot(width=0.3, alpha=0.7) +
    scale_fill_viridis(discrete = TRUE) +
    ylim(ylim)+
    theme_light() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("") +
    ylab(ylab) + 
    xlab(xlab)
  # coord_flip()
  
}


### "data": N_trial*N_sim_setting
plot_pointrange = function(data, ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE))){  
  library(ggplot2)
  library(dplyr)
  library(viridis)
  
  ylim = ylim
  
  data = reshape2::melt(data, id.vars=NULL, variable.name="name")
  method = rep("1",time=dim(data)[1])
  data %>%
    # left_join(method) %>%
    # mutate(myaxis = paste0(name)) %>%
    ggplot( aes(x=name, y=value, col=method, group=method)) +
    stat_summary(geom = "line",fun=median)+
    stat_summary(geom = "point",fun=median)+
    stat_summary(geom = "errorbar",width=0.1,
                 fun.min = function(z) { quantile(z,0.25) },
                 fun.max = function(z) { quantile(z,0.75) }
    )+
    coord_cartesian(ylim=ylim)+
    theme_light() +
    theme(
      axis.title = element_blank(),
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("") 
  # coord_flip()
  
}


plot_ARI_compr = function(ARI_our, ARI_ppsbm, ylim=c(0,1), continuous=FALSE, reverse=FALSE, n.breaks=10){
  data_tmp=rbind(reshape2::melt(ARI_our), reshape2::melt(ARI_ppsbm))
  method = rep(c("our method","ppsbm"),times=c(nrow(reshape2::melt(ARI_our)), nrow(reshape2::melt(ARI_ppsbm))))
  data_tmp = cbind(data_tmp, method)
  colnames(data_tmp) = c("V","ARI","method")
  
  if(continuous){
    data_tmp$V <- as.numeric(as.character(data_tmp$V))
    width = (max(data_tmp$V)-min(data_tmp$V))/50
  }
  else{
    width = 0.1
  }
  
  # grouped boxplot
  g<- ggplot(data_tmp, aes(x=V, y=ARI, group=method, col=method)) + 
    # geom_boxplot()+
    stat_summary(geom = "line",fun=median)+
    stat_summary(geom = "point",fun=median)+
    stat_summary(geom = "errorbar",width=width,
                 fun.min = function(z) { quantile(z,0.25) },
                 fun.max = function(z) { quantile(z,0.75) }
    )+
    coord_cartesian(ylim = ylim)+
    theme_light() +
    theme(
      axis.title = element_blank(),
      # legend.position="none",
      plot.title = element_text(size=11)
    ) #+
    # scale_x_discrete(guide = guide_axis(n.dodge=2, check.overlap = T))
  
  if(reverse&continuous)
    g + scale_x_reverse(breaks=unique(data_tmp$V))
  else if(continuous)
    g+scale_x_continuous(breaks=unique(data_tmp$V))
  else
    g
}




# Extract certain measurement from a file_list ---------------------------------


extract_measurement = function(file_list, measurement){
  measurement_df = c()
  for (file in file_list) {
    load(file)
    tmp_df = c()
    for (i in 1:length(results1)) {
      result_setting_i = results1[[i]]
      # browser()
      tmp = sapply(result_setting_i, function(one_trial) tryCatch(one_trial[[measurement]], 
                                                                  error=function(x)NA))
      tmp_df = cbind(tmp_df, tmp)
    }
    tmp_df = as.data.frame(tmp_df)
    measurement_df = rbind(measurement_df, tmp_df)
  }
  colnames(measurement_df) = names(results1)
  return(measurement_df)
}








