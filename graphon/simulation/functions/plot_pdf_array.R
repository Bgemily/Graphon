
plot_pdf_array = function(pdf_array_list, pdf_true_array = NULL, t_vec = seq(0, 50, 0.05))
{
  if (!is.list(pdf_array_list))
    pdf_array_list = list(pdf_array_list)
    
  t_unit = t_vec[2] - t_vec[1]
  
  k1 = dim(pdf_array_list[[1]])[1]
  k2 = dim(pdf_array_list[[1]])[2]
  # par(mfrow = c(k1,k2))
  
  # if(!is.null(pdf_true_array) & k1==k2){ 
  #   if (k1!=dim(pdf_true_array)[1] | k2!=dim(pdf_true_array)[2])
  #     stop("size of pdf_array and pdf_true_array should match.")
  #   
  #   # find permutation for pdf_array
  #   permn_list = lapply(pdf_array_list, function(pdf_array)find_permn(pdf_array, pdf_true_array, t_unit = t_unit)$permn)
  #   pdf_array_list = mapply(function(pdf_array, permn)pdf_array[permn,permn,], pdf_array_list, permn_list, SIMPLIFY = FALSE)
  # }
  # else if (k1==k2){
  #   permn_list = lapply(pdf_array_list, function(pdf_array)find_permn(pdf_array, pdf_true_array, t_unit = t_unit)$permn)
  #   pdf_array_list = mapply(function(pdf_array, permn)pdf_array[permn,permn,], pdf_array_list, permn_list, SIMPLIFY = FALSE)
  # }
  # 
  
  big.df = data.frame()
  for (q in 1:k1) {
    for (l in 1:k2) {
      
      pdf_array_mat = sapply(pdf_array_list, function(pdf_array)pdf_array[q,l,]) # length(t_vec)*N_subj
      tmp.df = data.frame(cbind(t=t_vec, true=pdf_true_array[q,l,], pdf_array_mat, mean=rowMeans(pdf_array_mat))) 
      tmp.df = reshape2::melt(tmp.df, id.vars = 't')
      tmp.df$col.group = sapply(as.vector(tmp.df$variable), switch, true="True", mean="Mean", "Estimate")

      tmp.df$clus.pair = paste0("(",q,",",l,")")
      
      big.df = rbind(big.df, tmp.df)
    }
  }
  
  library(ggplot2)
  ggplot(big.df, aes(t, value, group=variable, color=col.group, alpha=col.group, size=col.group, linetype=col.group)) +
    geom_line() + 
    scale_color_manual(values = c("Estimate"="black","True"="red", "Mean"="black"), name=NULL) + 
    scale_alpha_manual(values=c("Estimate"=0.2,"True"=1, "Mean"=1), name=NULL) +
    scale_size_manual(values=c("Estimate"=0.2,"True"=0.7, "Mean"=0.4), name=NULL)+
    scale_linetype_manual(values = c("Estimate"=1,"True"=1, "Mean"=2)) + 
    ylim(c(0, 0.4)) +
    facet_wrap(~clus.pair) +
    xlab("Time") + ylab(NULL) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(strip.background = element_blank(),  strip.text.x = element_blank()) + 
    geom_text(data=big.df[!duplicated(big.df$clus.pair),], x=Inf, y=Inf, aes(label=clus.pair), size=3, color='black',
              vjust = "inward", hjust = "inward", 
             show.legend = FALSE)
  

}
