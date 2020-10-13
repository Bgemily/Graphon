visual_realdata = function(result, path, suffix="_L.pdf")
{
  tmp = result
  
  # plot estimated f_qk's
  adjs_edge_time_mat = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                                          n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                                          zlim = NULL,reorder = FALSE, showplot = FALSE)
  center_pdf_array = get_center_pdf_array(edge_time_mat = adjs_edge_time_mat, clusters = tmp$clus_result$clusters, 
                                          n0_vec = tmp$clus_result$n0_vec, n0_mat = 0*tmp$clus_result$n0_mat, 
                                          t_vec = tmp$network$t_vec, bw = bw)
  pdf(paste0("./plots/",'conn_patt_',path,suffix),width = 4,height = 4)
  g <- plot_pdf_array(center_pdf_array,t_vec = tmp$network$t_vec,y_lim=c(0,0.055))
  print(g)
  dev.off()
  
  # plot edge time matrix heatmap (with or w/out activation time; with or w/out reordering)
  pdf(paste0("./plots/",'edge_time_denoise_reorder_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = TRUE)
  dev.off()
  
  pdf(paste0("./plots/",'edge_time_noise_reorder_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = TRUE)
  dev.off()
  
  pdf(paste0("./plots/",'edge_time_denoise_no_order_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = FALSE)
  dev.off()
  
  pdf(paste0("./plots/",'edge_time_noise_no_order_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = FALSE)
  dev.off()
  
  
  ### plot heatmap of activation time
  pdf(paste0("./plots/",'max_active_time_heatmap_',path,suffix),width = 4,height = 4)
  activ_time_mat = tmp$clus_result$n0_mat * tmp$network$t_vec[2]
  plot_edge_time_mat(edge_time_mat = activ_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = TRUE)
  dev.off()
  
  
  ### plot cluster membership matrix
  pdf(paste0("./plots/",'clus_mem_no_order_',path,suffix),width = 4,height = 4)
  image(as.matrix(clus2mem(tmp$clus_result$clusters)), col=gray.colors(3,start=1,end=0))
  dev.off()
  
  pdf(paste0("./plots/",'clus_mem_reorder_',path,suffix),width = 4,height = 4)
  image(as.matrix(sort(clus2mem(tmp$clus_result$clusters))), col=gray.colors(3,start=1,end=0))
  dev.off()
  
  mem_mat = dummies::dummy(clus2mem(tmp$clus_result$clusters))
  pdf(paste0("./plots/",'clus_mem_mat_',path,suffix),width = 4,height = 4)
  image(t(mem_mat), col=gray.colors(2,start = 1,end=0))
  dev.off()
  
  
  ### plot distribution of activation time
  data = data.frame(value=tmp$clus_result$n0_vec*tmp$network$t_vec[2], 
                    type=clus2mem(tmp$clus_result$clusters[]))
  pdf(paste0("./plots/",'active_time_hist_',path,suffix),width = 4,height = 4)
  g<- ggplot(data, aes(x=value,fill=as.factor(type)))+
    geom_histogram( aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins=5 ) +
    scale_fill_manual(values=2:4) +
    theme_bw() +
    labs(fill="")+
    facet_wrap(~type)
  print(g)
  dev.off()
  
}


visual_realdata_2 = function(result, path, suffix="_L.pdf", edge.time, locs, 
                             reduced.dFF, window_length, window_step, cor.full.ave, member.ship){
  tmp=result
  
  ### plot network development animation (or a snapshot)
  edge.time = tmp$network$edge_time_mat
  edge.time = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                                 n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                                 zlim = NULL,reorder = FALSE, showplot = FALSE)
  
  max_conn_time = max(edge.time[which(edge.time<Inf)])
  mins = lapply(seq(0,max_conn_time,5), function(t)c(0,t))
  
  # plot_network(locs = cbind(locs[tmp$id,2], -locs[tmp$id,1]), edge.time = edge.time, 
  #              output = paste0("animation_denoise",'_',path,suffix, ".gif"),
  #              window_list = list(mins), asp=2, save_plots = T, delay=20, 
  #              cols = t(col2rgb(1+member.ship[tmp$id]*0)), alpha=100)
  
  
  plot_network(locs = locs[,], edge.time = edge.time[,], 
               output = paste0("spatial_location_",path,suffix), vertex.size = 3,
               window_list = list(0), asp=0.3, save_plots = T, delay=20, 
               cols = t(col2rgb(1+member.ship[])))
  
  ### plot heatmap for neural activity traces
  n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
  interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
  interval.list = interval.list + (1:n.intervals-1)*window_step
  ave_dFF = matrix(nrow=dim(reduced.dFF)[1], ncol=n.intervals)
  for(i in 1:n.intervals){
    ave_dFF[,i] = rowMeans(reduced.dFF[,interval.list[i,]]);
  }
  
  pdf(paste0("./plots/",'activity_heatmap_',path,suffix),width = 4,height = 4)
  fields::image.plot(t(ave_dFF[tmp$id[unlist(tmp$clus_result$clusters)],]),zlim=c(0,0.05))
  dev.off()
  
  
  ### plot correlation curves
  pdf(paste0("./plots/",'corr_curves_',path,suffix),width = 4,height = 4)
  plot( (1:dim(cor.full.ave)[1])*window_step/240, 
       cor.full.ave[,sample(tmp$id[tmp$clus_result$clusters[[1]]],size = 1),
                    sample(tmp$id[tmp$clus_result$clusters[[1]]],size = 1) ],
       type='l',col=1,xlim=c(0,280),ylim=c(-0.2,1),ylab='',xlab='')
  abline(h=0.6,col=4,lty=3)
  for (q in 1:3) {
    for (k in q:3) {
      for (. in 1:5) {
        lines(cor.full.ave[,sample(tmp$id[tmp$clus_result$clusters[[q]]],size = 1),
                           sample(tmp$id[tmp$clus_result$clusters[[k]]],size = 1) ],
              type='l',col=ifelse(k==3&q!=2,1,2),xlim=c(0,280),ylim=c(-0.2,1),ylab='',xlab='')
      }
    }
  }
  dev.off()
  
  
}










