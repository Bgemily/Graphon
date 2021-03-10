
rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)



# Compare N_subj ----------------------------------------------------------
library("tidyverse")

path_vec = rep(0,3)
path_vec[1] = "../Results/Rdata/20210307_224827"
path_vec[2] = "../Results/Rdata/20210307_200804"
path_vec[3] = "../Results/Rdata/20210307_223906"

file_list_mat = sapply(path_vec, function(path)list.files(path=path, full.names = T))

Lambda_mean_sq_err_longlong = data.frame()
ARI_mean_longlong = data.frame()
for (j in 1:3) {
  file_list = file_list_mat[,j]
  Lambda_mean_sq_err = extract_measurement(file_list = file_list, measurement = "Lambda_mean_sq_err")
  Lambda_mean_sq_err_long = Lambda_mean_sq_err %>% 
                                pivot_longer(cols = starts_with("N_subj"), 
                                             names_to = "N_subj", 
                                             values_to = "Lambda_mse") %>%
                                mutate(Setup = j, N_subj=recode(N_subj, N_subj_1=1L, N_subj_5=5L, N_subj_10=10L))
  Lambda_mean_sq_err_longlong = rbind(Lambda_mean_sq_err_longlong, Lambda_mean_sq_err_long)
  
  ARI_mean = extract_measurement(file_list = file_list, measurement = "ARI_mean")
  ARI_mean_long = ARI_mean %>% 
                    pivot_longer(cols = starts_with("N_subj"), 
                                  names_to = "N_subj", 
                                  values_to = "ARI") %>%
                    mutate(Setup=j, N_subj=recode(N_subj, N_subj_1=1L, N_subj_5=5L, N_subj_10=10L)) 
  ARI_mean_longlong = rbind(ARI_mean_longlong, ARI_mean_long)
}

# pdf(file="../Results/Plots/Lambda_vs_N_subj.pdf", width = 4, height = 4)
Lambda_mean_sq_err_longlong %>%
  mutate(Setup=recode(Setup, `1`="Ideal case",
                      `2`="Smaller beta", 
                      `3`="Smaller beta & p")) %>%
  mutate(Setup=as.factor(Setup), N_subj=as.factor(N_subj)) %>%
  ggplot(mapping = aes(x=Setup, y=log(Lambda_mse), fill=N_subj)) +
  geom_boxplot() 
# dev.off()


# pdf(file="../Results/Plots/ARI_vs_N_subj.pdf", width = 4, height = 4)
ARI_mean_longlong %>%
  mutate(Setup=recode(Setup, `1`="Ideal case",
                      `2`="Smaller beta", 
                      `3`="Smaller beta & p")) %>%
  mutate(Setup=as.factor(Setup), N_subj=as.factor(N_subj)) %>%
  ggplot(mapping = aes(x=Setup, y=1-ARI, fill=N_subj)) +
  geom_boxplot() 
# dev.off()

# ARI ---------------------------------------------------------------------


ARI_ppsbm_0 = ARI_ppsbm_40 = c()
path_tmp = 'sim_results_clus_size/'
file_list = list.files(path=path_tmp, pattern='*ppsbm*',full.names = T)
for (file in file_list) {
  load(file)
  ARI = c()
  clus_size_mat = matrix(30,nrow=length(results1),ncol=3)
  clus_size_mat_tmp = t(sapply(clus_size_vec_list, function(x)x))
  clus_size_mat[(nrow(clus_size_mat)-nrow(clus_size_mat_tmp)+1):nrow(clus_size_mat),] = clus_size_mat_tmp
  for (i in 1:length(results1)) {
    result = results1[[i]]
    memb_true_vec = c(rep(1,clus_size_mat[i,1]), rep(2,clus_size_mat[i,2]),rep(3,clus_size_mat[i,3]))
    # tmp = sapply(result, function(tmp)tryCatch(get_one_ARI(memb_est_vec = clus2mem(tmp$clus_result$clusters),
    #                                                        memb_true_vec = memb_true_vec), error=function(x)Inf) )
    # tmp = sapply(result, function(tmp)tryCatch(tmp$clus_result$risk, error=function(x)0))
    tmp = sapply(result, function(tmp)tryCatch(tmp$clus_result$ARI, error=function(x)0))
    ARI = cbind(ARI, tmp)
  }
  ARI = as.data.frame(ARI)
  if(tau_std==0){
    ARI_ppsbm_0 = rbind(ARI_ppsbm_0,ARI)
  }
  else if(tau_std==40){
    ARI_ppsbm_40 = rbind(ARI_ppsbm_40, ARI)
  }
}
colnames(ARI_ppsbm_0) = colnames(ARI_ppsbm_40) = c(tau_max_vec,conn_prob_vec, beta_vec, alpha_vec,
                                                   sapply(clus_size_vec_list, function(x)x[3]))
# colnames(ARI_ppsbm_0) = colnames(ARI_ppsbm_40) = sapply(clus_size_vec_list, function(x)x[3])
plot_pointrange(data=ARI_ppsbm_40[,],ylim=c(0,1))
save(list=c("ARI_ppsbm_0","conn_prob_std","beta_std","alpha_std"),file=paste0(path_tmp,"ARI_psbm_0.rdata"))
save(list=c("ARI_ppsbm_40","conn_prob_std","beta_std","alpha_std"),file=paste0(path_tmp,"ARI_psbm_40.rdata"))
rm(path_tmp)




ARI_our=ARI_our_40; ARI_ppsbm=ARI_ppsbm_40
plot_ARI_compr(ARI_our[], ARI_ppsbm[], ylim=c(0,1))

plot_ARI_compr(ARI_our[,-7], ARI_ppsbm[,-7], ylim=c(0,1))

plot_ARI_compr(ARI_our[,c(1:5,9,12)], ARI_ppsbm[,c(1:5,9,12)], ylim=c(0,1))

plot_ARI_compr(ARI_our[,c(1:9)], ARI_ppsbm[,c(1:9)], ylim=c(0,1), 
               continuous = T, reverse=F)
plot_ARI_compr(ARI_our[,c(10:17)], ARI_ppsbm[,c(10:17)])
plot_ARI_compr(ARI_our[,c(18:23)], ARI_ppsbm[,c(18:23)], ylim=c(0,1))
plot_ARI_compr(ARI_our[,c(24:30)], ARI_ppsbm[,c(24:30)], ylim=c(0,1), continuous = T)
plot_ARI_compr(ARI_our[,c(31:41)], ARI_ppsbm[,c(31:41)], ylim=c(0,1))
plot_ARI_compr(ARI_our[,c(31,36:38)], ARI_ppsbm[,c(31,36:38)], ylim=c(0,1))



###  computing time
cluster_time = sapply(results1, function(result)sapply(result, function(x)as.double.difftime(x$clus_result$cluster_time, units="secs")))
align_time = sapply(results1[], function(result)sapply(result, function(x)as.double.difftime(x$clus_result$align_time, units = "secs")))
N_iteration = sapply(results1[],function(result)sapply(result,function(x)tryCatch(length(x$clus_result$clusters_history), error=function(x)11)))
cluster_time = colMeans(cluster_time)
align_time = colMeans(align_time)
running_time = cluster_time + align_time
N_iteration = colMeans(N_iteration)

plot(-log10(conv_threshold_vec[c(1,2,3,6,7,8)]),running_time[c(1,2,3,6,7,8)],type='b')
plot_pointrange(data=ARI[,])


### l2 distance (risk)
risk = c()
for (result in results1) {
  true_pdf_array = fun2pdfarray(true_pdf_fun_list = result[[1]]$network$true_pdf_fun_list,
               tau_mat = result[[1]]$network$tau_mat*0, 
               membership_true = result[[1]]$network$membership_true,
               t_vec = result[[1]]$network$t_vec)
  t_unit = result[[1]]$network$t_vec[2]
  conn_prob = sum(true_pdf_array[1,1,]*t_unit)
  N_clus = length(result[[1]]$clus_result$clusters)
  
  risk_vec = sapply(result, function(x)sqrt(sum((true_pdf_array-x$clus_result$center_pdf_array)^2)*t_unit/N_clus^2)/conn_prob)
  risk = cbind(risk, risk_vec)
}
risk = as.data.frame(risk)
colnames(risk) = names(results1)
colnames(risk) = c(tau_max_vec,conn_prob_vec, beta_vec)
plot_jitter_boxplot(data=risk[,], ylim = c(0,0.4))



plot_pointrange(risk_our[,c(1,2,4,6,7)], ylim=c(0,0.3))
plot_pointrange(risk_our[,c(8:11)], ylim=c(0,0.3))
plot_pointrange(risk_our[,c(13:17)], ylim=c(0,0.3))





# Plot estimated connecting patterns (simulation) -------------------------------------

results = results1$tau_max_20


pdf_true_array = fun2pdfarray(true_pdf_fun_list = results[[1]]$network$true_pdf_fun_list,
                               tau_mat = results[[1]]$network$tau_mat, 
                              membership_true = results[[1]]$network$membership_true,
                              t_vec = results[[1]]$network$t_vec)
# times pdf_true_array by connecting probabilities
# r = results[[2]]
# edge_time_mat = r$network$edge_time_mat
# clus_degree = get_clus_degree_mat(edge_time_mat = edge_time_mat, clusters = mem2clus(r$network$membership_true))
# for (q in 1:dim(pdf_true_array)[1]) {
#   for (l in 1:dim(pdf_true_array)[2]) {
#     pdf_true_array[q,l, ] = pdf_true_array[q,l, ] * clus_degree[q,l] / (sum(r$network$membership_true==q)*sum(r$network$membership_true==l)*(1-0.5*I(q==l)))
#   }
# }


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





# Test simulation setting's difficulty -------------------------------------

tmp = results1$`tau_max_20`[[1]]

main(case=2,SEED=9,N_clus=3,tau_max=c(40),conn_prob = 0.3,beta = 1.3, alpha=1,
     standardize = T, tau_struc = max, MaxIter = 10, total_time = 200)->tmp
tmp$clus_result$clusters_history
# tmp$clus_result$n0_vec*tmp$network$t_vec[2]
# plot_aggr_traces(tmp$network$edge_time_mat,color=tmp$network$membership_true,bw=2)
# plot_aggr_traces(tmp$network$edge_time_mat,color=clus2mem(tmp$clus_result$clusters),bw=2)
pdf_array = fun2pdfarray(true_pdf_fun_list = tmp$network$true_pdf_fun_list, tau_mat = tmp$network$tau_mat*0, 
                         membership_true = tmp$network$membership_true, t_vec = tmp$network$t_vec)
# plot_pdf_array(pdf_array,t_vec = tmp$network$t_vec,y_lim = c(0,0.03))
plot_pdf_array(pdf_array_list = tmp$clus_result$center_pdf_array, pdf_true_array = pdf_array,
               t_vec = tmp$network$t_vec,y_lim = c(0,0.03))
# plot_pdf_array(pdf_array_list = tmp$clus_result$center_cdf_array, pdf_true_array = NULL,t_vec = tmp$network$t_vec)


plot(tmp$network$tau_vec, tmp$clus_result$n0_vec*tmp$network$t_vec[2]); abline(a=0,b=1,col=2)
earliest_edge_time = apply(tmp$network$edge_time_mat, 1, function(row)min(row[which(row>1)]))
plot(tmp$network$tau_vec,earliest_edge_time); abline(a=0,b=1,col=2)


apply_ppsbm(case=2,SEED=81, Qmin=3, 
            tau_max=40,conn_prob = 0.3, beta=1.3, alpha=1,
            tau_struc = max,total_time = 200)->tmp2
tmp2$clus_result$clusters

plot_aggr_traces(tmp2$network$edge_time_mat,color=clus2mem(tmp2$clus_result$clusters),bw=10,xlim=c(0,340),ylim = c(0,0.04))
par(mfrow=c(2,3))
for (ql in 1:nrow(tmp2$clus_result$logintensities.ql)) {
  plot(exp(tmp2$clus_result$logintensities.ql[ql,]),type='l')
}
par(mfrow=c(1,1))






# plot connecting patterns obtained by ppsbm ------------------------------

tmp = results1$conn_prob_0.2[[sample(100,1)]]
intensity_tmp = tmp$clus_result$logintensities.ql
intensity_tmp[intensity_tmp==0] = -Inf
intensity_tmp = exp(intensity_tmp)

z = rbind(rep(c(1,0,0),each=30),
          rep(c(0,1,0),each=30),
          rep(c(0,0,1),each=30))
mem_tmp = apply(tmp$clus_result$tau, 2, which.max)
intensity_tmp = ppsbm::sortIntensities(intensity_tmp,z,tmp$clus_result$tau, directed = F)

ind_tmp = matrix(c(1,2,3,2,4,5,3,5,6),nrow=3,byrow = T)
for (q in 1:3) {
  for (l in 1:3) {
    hist_tmp = density(tmp$network$edge_time_mat[which(mem_tmp==q),which(mem_tmp==l)])
    prob_tmp = sum(tmp$network$edge_time_mat[which(mem_tmp==q),which(mem_tmp==l)]<Inf)/(sum(mem_tmp==q)*sum(mem_tmp==l))
    hist_tmp$y = hist_tmp$y * prob_tmp
    i = ind_tmp[q,l]
    plot(seq(0,200,length.out = dim(intensity_tmp)[2]),intensity_tmp[i,],type='l',col=2)
    # lines(hist_tmp,xlim = c(0,200))
    
  }
}

tmp$clus_result$clusters


# Real data result  -------------------------------------------------------


tmp=L_result
tmp=R_result
temp=c(1,2,3)
plot_pdf_array(tmp$clus_result$center_pdf_array[temp,temp,],t_vec = tmp$network$t_vec,y_lim=c(0,0.05))
tmp$clus_result$clusters
# tmp$clus_result$L_mat
tmp$clus_result$n0_vec[tmp$clus_result$clusters[[3]]]*tmp$network$t_vec[2]
edge_time_mat[tmp$id[tmp$clus_result$clusters[[3]]],
              tmp$id[tmp$clus_result$clusters[[3]]]]
tmp$clus_result$n0_mat[tmp$clus_result$clusters[[3]],tmp$clus_result$clusters[[2]]]*tmp$network$t_vec[2]


# plot estimated v_i vs patterning time
earliest_edge_time = apply(tmp$network$edge_time_mat, 1, function(row)min(row[which(row>1)]))
plot(tmp$clus_result$n0_vec*tmp$network$t_vec[2], earliest_edge_time); abline(a=0,b=1,col=2)


data_folder = "../processed_FunctionalData/"
path.list=list.files(data_folder);

# plot estimated connecting patterns (real data) -----
for(k in 1:length(path.list)){
  path = path.list[[k]]
  tmp = load(list.files(pattern = path))
  sapply(file.sources, source)
  tmp=L_result
  pdf(paste0("./plots/",path,"_Lside.pdf"),width = 4,height = 4)
  adjs_edge_time_mat = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                                          n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                                          zlim = NULL,reorder = FALSE, showplot = FALSE)
  center_pdf_array = get_center_pdf_array(edge_time_mat = adjs_edge_time_mat, clusters = tmp$clus_result$clusters, 
                                          n0_vec = tmp$clus_result$n0_vec, n0_mat = 0*tmp$clus_result$n0_mat, 
                                          t_vec = tmp$network$t_vec, bw = bw)
  plot_pdf_array(center_pdf_array,t_vec = tmp$network$t_vec,y_lim=c(0,0.05))
  dev.off()
  tmp=R_result
  pdf(paste0("./plots/",path,"_Rside.pdf"),width = 4,height = 4)
  adjs_edge_time_mat = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                                          n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                                          zlim = NULL,reorder = FALSE, showplot = FALSE)
  center_pdf_array = get_center_pdf_array(edge_time_mat = adjs_edge_time_mat, clusters = tmp$clus_result$clusters, 
                                          n0_vec = tmp$clus_result$n0_vec, n0_mat = 0*tmp$clus_result$n0_mat, 
                                          t_vec = tmp$network$t_vec, bw = bw)
  plot_pdf_array(center_pdf_array,t_vec = tmp$network$t_vec,y_lim=c(0,0.05))
  dev.off()
}


tmp=L_result_new
islet = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/islet.csv',sep='')))
islet = islet[,-1]

file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

adjs_edge_time_mat = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                               n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec,
                               zlim = NULL,reorder = FALSE, showplot = FALSE)
center_pdf_array = get_center_pdf_array(edge_time_mat = adjs_edge_time_mat, clusters = tmp$clus_result$clusters, 
                                        n0_vec = tmp$clus_result$n0_vec, n0_mat = 0*tmp$clus_result$n0_mat, 
                                        t_vec = tmp$network$t_vec, bw = bw)
center_pdf_array = center_pdf_array[c(2:4,1),c(2:4,1),]
plot_pdf_array(center_pdf_array[],t_vec = tmp$network$t_vec,y_lim=c(0,0.055))


##### Calculate islet+/-/NA composition of clusters -----

get_composition_array = function(file_list, Lside=T){
  array_tmp = array(dim = c(4,7,length(file_list)))
  for (i in 1:length(file_list)) {
    file = file_list[i]
    load(file)
    if(Lside){
      tmp=L_result_new
    }
    else{
      tmp=R_result_new
    }
    composition_matrix = matrix(0,nrow=length(tmp$clus_result$clusters), ncol=7)
    for (q in 1:length(tmp$clus_result$clusters)) {
      tmp_L = L_result_new; tmp_R = R_result_new
      islet_res = islet[c(tmp_L$id[tmp_L$clus_result$clusters[[q]]],tmp_R$id[tmp_R$clus_result$clusters[[q]]])]
      composition_matrix[q,2:4] = c(length(which(islet_res==1)), length(which(islet_res==0)), sum(is.na(islet_res))) 
      composition_matrix[q,1] = length(islet_res)
      composition_matrix[q,5:7] = composition_matrix[q,2:4]/composition_matrix[q,1]
    }
    colnames(composition_matrix) = c("Cluster Size", "MN", "Ventral IN", "CoLA/CoSA", "MN(%)", "Ventral IN(%)", "CoLA/CoSA(%)")
    rownames(composition_matrix) = c("Cluster 2", "Cluster 1", "Cluster 3", "Iso Nodes")[1:nrow(composition_matrix)]
    array_tmp[1:nrow(composition_matrix),,i]=composition_matrix
  }
  dimnames(array_tmp) <- list(  rownames(composition_matrix), colnames(composition_matrix) ,seq(0.1,0.7,0.1))
  return(array_tmp[c(2,1,3,4),,])
}

file_list = list.files(path='/Users/bgemily/Dropbox/DynamicSBM/Data/real_data_results/islet_20170202/',pattern = "*Nclus3",full.names = T)
array_tmp_1 = get_composition_array(file_list, Lside=T)

file_list = list.files(path='/Users/bgemily/Dropbox/DynamicSBM/Data/real_data_results/islet_20170216/',pattern = "*Nclus3",full.names = T)
array_tmp_2 = get_composition_array(file_list, Lside=F)


thres = 4

## composition of estimated clusters: scatter plot
plot(1,type='n',xlim=c(0,1)-0.5,ylim=c(0,1)-0.5,xlab="relative MN(%)",ylab="relative Ventral IN(%)", main="cluster composition")
abline(h=0,lty=2,col="gray"); abline(v=0,lty=2,col="gray")
overall_MN = sum(array_tmp_1[,2,thres])/sum(array_tmp_1[,1,thres])
overall_VenIN = sum(array_tmp_1[,3,thres])/sum(array_tmp_1[,1,thres])
for (q in 1:4) {
  for (thres in thres) {
    points(array_tmp_1[q,5,thres]-overall_MN, array_tmp_1[q,6,thres]-overall_VenIN, 
           col=switch(q,"red","green","blue","black"),
           cex=(array_tmp_1[q,1,thres]/10),
           pch=21)
  }
}
overall_MN = sum(array_tmp_2[,2,thres])/sum(array_tmp_2[,1,thres])
overall_VenIN = sum(array_tmp_2[,3,thres])/sum(array_tmp_2[,1,thres])
for (q in 1:4) {
  for (thres in thres) {
    points(array_tmp_2[q,5,thres]-overall_MN, array_tmp_2[q,6,thres]-overall_VenIN, 
           col=switch(q,"red","green","blue","black"),
           cex=(array_tmp_2[q,1,thres]/10),
           pch=24)
  }
}

# composition of estimated clusters: stacked barplot
array_tmp = apply(array_tmp_1,c(1,3),function(x)x[2:4]/sum(x[2:4]))[,,thres]
clus_size_perc = colSums(array_tmp_1[,2:4,thres])/sum(array_tmp_1[,1,thres])
array_tmp = cbind(array_tmp,Marginal=clus_size_perc)
# array_tmp = sweep(array_tmp,1,clus_size_perc)

data = reshape2::melt(array_tmp)
names(data) = c("Cell_type","Cluster","Percentage")
library(ggplot2)
ggplot(data, aes(fill=Cell_type, y=Percentage, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

array_tmp = apply(array_tmp_2,c(1,3),function(x)x[2:4]/sum(x[2:4]))[,,thres]
clus_size_perc = colSums(array_tmp_2[,2:4,thres])/sum(array_tmp_2[,1,thres])
array_tmp = cbind(array_tmp,Marginal=clus_size_perc)
# array_tmp = sweep(array_tmp,1,clus_size_perc)

data = reshape2::melt(array_tmp)
names(data) = c("Cell_type","Cluster","Percentage")
library(ggplot2)
ggplot(data, aes(fill=Cell_type, y=Percentage, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")



## composition of cell types
array_tmp = apply(array_tmp_1,c(2,3),function(x)x/sum(x))[,2:4,thres]
clus_size_perc = array_tmp_1[,1,thres]/sum(array_tmp_1[,1,thres])
array_tmp = cbind(array_tmp,Marginal=clus_size_perc)
# array_tmp = sweep(array_tmp,1,clus_size_perc)

data = reshape2::melt(array_tmp)
names(data) = c("Cluster","cell_type","percentage")
library(ggplot2)
ggplot(data, aes(fill=Cluster, y=percentage, x=cell_type)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Cell type") + ylab("Percentage")

array_tmp = apply(array_tmp_2,c(2,3),function(x)x/sum(x))[,2:4,thres]
clus_size_perc = array_tmp_2[,1,thres]/sum(array_tmp_2[,1,thres])
array_tmp = cbind(array_tmp,Marginal=clus_size_perc)
# array_tmp = sweep(array_tmp,1,clus_size_perc)

data = reshape2::melt(array_tmp)
names(data) = c("Cluster","cell_type","value")
library(ggplot2)
ggplot(data, aes(fill=Cluster, y=value, x=cell_type)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Cell type") + ylab("Percentage")


array_tmp = apply(array_tmp_1+array_tmp_2,c(2,3),function(x)x/sum(x))[,2:4,thres]
colSums((array_tmp_1+array_tmp_2)[,2:4,thres])
clus_size_perc = (array_tmp_1+array_tmp_2)[,1,thres]/sum((array_tmp_1+array_tmp_2)[,1,thres])
array_tmp = cbind(array_tmp,Marginal=clus_size_perc)
# array_tmp = sweep(array_tmp,1,clus_size_perc)

data = reshape2::melt(array_tmp)
names(data) = c("Cluster","cell_type","value")
library(ggplot2)
ggplot(data, aes(fill=Cluster, y=value, x=cell_type)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Cell type") + ylab("Percentage")+
  scale_x_discrete(labels=c("MN" = "MN(n=73)", "Ventral IN" = "Ventral IN(n=39)",
                            "CoLA/CoSA" = "CoLA/CoSA(n=19)"))


# plot edge time matrix heatmap (with or w/out activation time; with or w/out reordering) -----
tmp=L_result
tmp=R_result
plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = TRUE)
plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = TRUE)
plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = FALSE)
plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = FALSE)


### plot heatmap of activation time -----
tmp = L_result

activ_time_mat = tmp$clus_result$n0_mat * tmp$network$t_vec[2]

plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat-activ_time_mat+30, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = TRUE)
plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = TRUE)

plot_edge_time_mat(edge_time_mat = activ_time_mat, clusters = tmp$clus_result$clusters,
                   n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                   zlim = c(0,300),reorder = TRUE)



### plot cluster membership matrix -----
tmp = L_result
image(as.matrix(clus2mem(tmp$clus_result$clusters)), col=gray.colors(3,start=1,end=0))
image(as.matrix(sort(clus2mem(tmp$clus_result$clusters))), col=gray.colors(3,start=1,end=0))
mem_mat = dummies::dummy(clus2mem(tmp$clus_result$clusters))
image(t(mem_mat), col=gray.colors(2,start = 1,end=0))

### plot network development animation (or a snapshot) ------
data_folder = "../processed_FunctionalData/"
path.list=list.files(data_folder);

for (k in 1:length(path.list)) {
  path=path.list[[k]]
  
  tmp = load(list.files(pattern = path))
  
  avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
  avai.inds=avai.inds[,-1];
  
  locs.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs=locs.all[avai.inds,]
  
  member.ship = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/MembShip.csv',sep='')))
  member.ship=member.ship[,-1];
  
  mnx.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/mnx.csv',sep='')))
  mnx.all = mnx.all[,-1]
  mnx = mnx.all[avai.inds]
  
  tmp=R_result
  edge.time = tmp$network$edge_time_mat
  edge.time = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                     zlim = NULL,reorder = FALSE, showplot = FALSE)
  
  max_conn_time = max(edge.time[which(edge.time<Inf)])
  mins = lapply(seq(70,130,5), function(t)c(0,t))
  
  member.ship=membership
  member.ship[membership==1]=2
  member.ship[membership==2]=1
  palette("default")
  plot_network_animation(locs = cbind(locs[tmp$id,2], -locs[tmp$id,1]), edge.time = edge.time, 
               output = "animation_denoise.gif",
               window_list = mins, asp=2, save_plots = F, delay=100, 
               cols = t(col2rgb(1+member.ship[tmp$id])), alpha=100)
  
  
  edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
  edge.time=edge.time[,-1]
  plot_network(locs = locs[,], edge.time = edge_time_mat[], output = "./plots/",
               filename = "Network2.pdf", vertex.size = 3,
               window_list = list(0), asp=0.5, save_plots = T, delay=20, 
               cols = t(col2rgb(1+member.ship[])))
  
  
  order.temp = order(locs[,1], decreasing = TRUE)
  order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])
  
}


### plot distribution of activation time ------
tmp = R_result_new
# path=path.list[[k]]
# member.ship = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/MembShip.csv',sep='')))
# member.ship=member.ship[,-1];

data = data.frame(value=tmp$clus_result$n0_vec*tmp$network$t_vec[2], type=clus2mem(tmp$clus_result$clusters[]))
# data = data.frame(value=c(L_result$clus_result$n0_vec,R_result$clus_result$n0_vec)*R_result$network$t_vec[2], 
#                   type=c(clus2mem(L_result$clus_result$clusters[]),
#                          clus2mem(R_result$clus_result$clusters[])))

data = data[which(data$type<4),]
data$type[which(data$type==1)] = 0
data$type[which(data$type==2)] = 1
data$type[which(data$type==0)] = 2
clus_size_vec = sapply(tmp$clus_result$clusters[c(2,1,3,4)], length)
data$clus_size = clus_size_vec[data$type]


library(ggplot2)
ggplot(data, aes(x=value,color=as.factor(type), size=as.factor(type)))+
  # geom_histogram( aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins=7 ) +
  geom_density( aes(), position = 'identity', bw=20) +
  scale_size_manual(values=clus_size_vec/10) +
  theme_bw() +
  theme(legend.position = "none")+
  xlab("Time(min)")+
  xlim(0,250)+
  ylab("Density")
  # labs(fill="")#+
  # facet_wrap(~type)






### plot heatmap for neural activity traces ------
data_folder = "../processed_FunctionalData/"
path.list=list.files(data_folder);

for (k in 1:length(path.list)) {
  path=path.list[[k]]
  
  avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
  avai.inds=avai.inds[,-1];
  
  locs.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs=locs.all[avai.inds,]
  
  member.ship = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/MembShip.csv',sep='')))
  member.ship=member.ship[,-1];
  
  mnx.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/mnx.csv',sep='')))
  mnx.all = mnx.all[,-1]
  mnx = mnx.all[avai.inds]
  
  library(data.table)
  dat.dFF=as.matrix(fread(paste('../processed_FunctionalData/',path,'/dFF.csv',sep='')))
  dat.dFF=dat.dFF[,-1]
  # Simply calculate the covariance among dFF within the time periods
  # Obtain the list of missing traces 
  na.inds= is.na(dat.dFF[,1]);
  reduced.dFF=dat.dFF[!na.inds,];
  
  window_length = 240 # = 1min at 4 Hz
  window_step = 240 # = 1min at 4Hz
  
  n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
  interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
  interval.list = interval.list + (1:n.intervals-1)*window_step
  
  ave_dFF = matrix(nrow=dim(reduced.dFF)[1], ncol=n.intervals)
  for(i in 1:n.intervals){
    ave_dFF[,i] = rowMeans(reduced.dFF[,interval.list[i,]]);
  }
  
  tmp = load(list.files(pattern = path))
  tmp = R_result
  fields::image.plot(t(ave_dFF[tmp$id[unlist(tmp$clus_result$clusters[1])],]),zlim=c(0,0.05), legend.shrink = 0.1,xaxt='n',yaxt='n')
  fields::image.plot(t(ave_dFF[sample(tmp$id[order(locs[tmp$id, 1], decreasing = TRUE)],10) ,]),zlim=c(0,0.05), xaxt='n',yaxt='n')
  image(t(ave_dFF[sample(tmp$id[order(locs[tmp$id, 1], decreasing = TRUE)], 10) ,]),zlim=c(0,0.05), col=fields::tim.colors(300), xaxt='n',yaxt='n')
  
}

fields::image.plot(t(ave_dFF[which(islet==1),]),zlim=c(0,0.05), legend.shrink = .5,xaxt='n',yaxt='n')


### plot edge times as lollipop plots -----

i=8
edge_time_vec = tmp$network$edge_time_mat[tmp$clus_result$clusters[[2]][i],
                                          unlist(tmp$clus_result$clusters)]
active_time = tmp$clus_result$n0_mat[tmp$clus_result$clusters[[2]][i],
                                     unlist(tmp$clus_result$clusters)]*tmp$network$t_vec[2]
active_time = active_time[which(edge_time_vec<Inf & edge_time_vec>0.1)]
edge_time_vec = edge_time_vec[which(edge_time_vec<Inf & edge_time_vec>0.1)]
data = data.frame(
  # x = as.integer(edge_time_vec-active_time),
  x = as.integer(edge_time_vec),
  y = 1
)
ggplot(data, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color=fields::tim.colors(300)[data$x], size=1) +
  # geom_point( color=fields::tim.colors(300)[data$x], size=4) +
  theme_void() +
  xlim(0,300)
ggsave("~/Dropbox/DynamicSBM/Data/plots/edge_time_lollipop_4.pdf",width = 10,height = 0.5)



### plot all correlation curves (with some representatives with solid colors) ------

id_tmp = which(locs[,2]<0)
plot(1, type="n", xlab="Time (min)", ylab="Correlation", xlim=c(0, 280), ylim=c(-0.1, 1))

for (i in 1:(length(id_tmp)-1)) {
  for (j in (i+1):length(id_tmp)) {
    cor_tmp = cor.full.ave[,id_tmp[i],id_tmp[j]]
    lines(cor_tmp, type='l',col=rgb(0,0,0,0.05),lwd=0.5) 
  }
}

for (i in 1:6) {
  cor_tmp = cor.full.ave[, sample(id_tmp,1), sample(id_tmp,1)]
  lines(cor_tmp, type='l',col=1,lwd=2) 
}


### plot maximum correlation heatmap -----

max_cor_mat = apply(cor.full.ave, c(2,3), max)
tmp = L_result_new
id_tmp = tmp$id[unlist(tmp$clus_result$clusters[])]
max_cor_mat = max_cor_mat[id_tmp,id_tmp]
fields::image.plot(max_cor_mat,zlim=c(0,1),col = fields::tim.colors(100),xaxt='n',yaxt='n')




### plot proportion vs time since passing threshold -----

pass_proportion = function(cor.full.ave, rho){
  exceed_time = apply(cor.full.ave, c(2,3), function(x)min(which(x>rho)))
  # drop_time = apply(cor.full.ave, c(2,3), function(x)ifelse(max(x)>rho & rev(x)[1]<rho,
  #                                                           max(which(x>rho)),
  #                                                           0))
  drop_time = exceed_time
  for (i in 1:dim(cor.full.ave)[2]) {
    for (j in 1:dim(cor.full.ave)[3]) {
      cor_tmp = cor.full.ave[,i,j]
      drop_time[i,j] = 1+min(which(cor_tmp[1:(length(cor_tmp)-1)]>=rho & cor_tmp[2:(length(cor_tmp))]<rho))
    }
  }
  
  id_tmp = which(exceed_time<Inf)
  drop_time[id_tmp] = drop_time[id_tmp] - exceed_time[id_tmp]
  # drop_time[drop_time<0] = Inf
  total = sum(exceed_time<Inf)
  
  prop_vec = rep(0,dim(cor.full.ave)[1])
  for (t in 1:dim(cor.full.ave)[1]) {
    prop_vec[t] = sum( drop_time<=t )
    prop_vec[t] = prop_vec[t] / total
  }
  
  id_2 = which(drop_time<75,arr.ind = T); # nodes that pass threshold but quickly drop back
  id_3 = which(drop_time>75&drop_time<200,arr.ind = T) # pass, and drop back after a while
  id_4 = which(exceed_time<Inf&drop_time==Inf,arr.ind = T) # pass, and not drop back
  id_1 = which(exceed_time==Inf,arr.ind = T) # never pass
  
  return(list(prop_vec=prop_vec,id_1=id_1,id_2=id_2, id_3=id_3,id_4=id_4))
}


### plot randomly selected correlation curves, colored by (1) never pass, (2) pass but quickly drop back,
### (3) pass, and drop back after a while, (4) pass, not drop back.
rho = 0.6
id_vec_tmp = which(locs[,2]<0)
res = pass_proportion(cor.full.ave = cor.full.ave[,id_vec_tmp,id_vec_tmp], rho=rho)
prop_vec = res$prop_vec
id_1 = res$id_1; id_2 = res$id_2; id_3 = res$id_3; id_4 = res$id_4
id_list = list(id_1,id_2,id_3,id_4)

plot(1, type="n", xlab="Time", ylab="Correlation", xlim = c(0,300),ylim=c(0,1))
for (. in 1:20) {
  group = sample(4,1)
  # group = 4
  id_tmp = id_list[[group]][sample(nrow(id_list[[group]]),1), ]
  lines(cor.full.ave[,id_vec_tmp,id_vec_tmp][,id_tmp[1],id_tmp[2]],col=group)
  # lines(cor.full.ave[,sample(id_5,1),sample(id_5,1)],col=5)
}
abline(h=rho,col=4,lty=2)
legend(x="topright",y=0.9, inset = c(-0.3,0), legend=c("Never pass","Drop quickly","Stay before drop","Never drop"),
       col=1:4, lty=1, cex=0.8)
#### ### ### ###

tmp = L_result

rho_vec = seq(0,0.9,0.1)
prop_mat = matrix(nrow=length(rho_vec),ncol=length(prop_vec))
for (i in 1:length(rho_vec)) {
  rho = rho_vec[i]
  prop_mat[i,] = pass_proportion(cor.full.ave = cor.full.ave[,tmp$id,tmp$id], rho=rho)$prop_vec
}

# plot proportion vs time since passing threshold for various rho
plot(1, type="n", xlab="Time since passing threshold", ylab="Proportion", xlim=c(0, 30), ylim=c(0, 1))
for (i in 1:length(rho_vec)) {
  lines(prop_mat[i,],col=i)
}
legend(x="topright",y=0.9, legend=rho_vec,
       col=1:(length(rho_vec)), lty=1, cex=0.8)

# plot proportion vs rho for different droping time
plot(1, type="n", xlab="rho", ylab="Proportion", xlim=c(0, 1), ylim=c(0, 1))
for (t in seq(1,20,1)) {
  lines(rho_vec,prop_mat[,t],col=t)
}
legend(x="topright",y=0.9, legend=seq(1,20,1),
       col=1:20, lty=1, cex=0.8)


# plot (proportion of quickly drop back nodes) vs rho for various max(cor) -------

quick_drop_proportion = function(cor.full.ave,rho,cor_max){
  max_cor_mat = apply(cor.full.ave, c(2,3), max)
  total = sum(max_cor_mat>cor_max & max_cor_mat<cor_max+0.1)
  N_quick_drop = 0
  for (i in 1:dim(cor.full.ave)[2]) {
    for (j in 1:dim(cor.full.ave)[3]) {
      if(max_cor_mat[i,j]>cor_max & max_cor_mat[i,j]<(cor_max+0.1)){
        exceed_time = min(which(cor.full.ave[,i,j]>rho))
        if(min(cor.full.ave[exceed_time:min(exceed_time+10,dim(cor.full.ave)[1]),i,j])<rho){
          N_quick_drop = N_quick_drop+1
        }
      }
    }
  }
  return(list(prop=N_quick_drop/total, total=total))
}

id_tmp = which(locs[,2]<0)
rho_vec = seq(0.1,0.9,0.1)
prop_mat = matrix(nrow=length(rho_vec),ncol=30)
total_vec = rho_vec
for (i in 1:length(rho_vec)) {
  cor_max = rho_vec[i]
  for (j in seq(30)) {
    rho = seq(0,cor_max,length.out=30)[j]
    res = quick_drop_proportion(cor.full.ave = cor.full.ave[,id_tmp,id_tmp], rho=rho, cor_max=cor_max)
    prop_mat[i,j] = res$prop
  }
  total_vec[i] = res$total
}

total_vec

max_cor_mat = apply(cor.full.ave, c(2,3), max)
hist(max_cor_mat[id_tmp,id_tmp], xlab="Maximum correlation", main = "")

par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=NA)
plot(1, type="n", xlab="Threshold", ylab="Proportion of quicly drop back correlations", xlim=c(0, 1), ylim=c(0, 1))
for (i in 1:length(rho_vec)) {
  cor_max = rho_vec[i]
  lines(seq(0,cor_max,length.out=30),prop_mat[i,],col=colorRamps::primary.colors(10)[i])
}
legend("topright", inset=c(-0.35,0), legend=rho_vec,
       col=colorRamps::primary.colors(10)[1:9], lty=1, cex=0.8)


autoimage::reset.par()





# Plot cluster memberships changes as rho grows ---------------------------

file_list = list.files(path = './real_data_results/func_20150410/', pattern = "*Nclus3_win240_rho",full.names = T)
membership_mat = c()
for (file in file_list) {
  load(file)
  membership[which(membership==0)] = 4
  membership_mat = rbind(membership_mat,membership[which(locs[,2]<0)])
}

rho_vec = seq(0.3,0.8,0.1)
clus_label_mat = matrix(c(3,1,2,4, 3,2,1,4, 2,3,1,4, 1,2,3,4, 2,1,3,4, 1,2,3,4),ncol=4,byrow=T)


# plot cluster membership migration using heatmap
N_nodes_mat_collapse = matrix(0,nrow=4,ncol=4)
for (from in 1:(length(file_list)-1)) {
  to = from+1
  N_nodes_mat = matrix(nrow=4,ncol=4)
  for (q in 1:4) {
    for (k in 1:4) {
      N_nodes_mat[q,k] = sum(membership_mat[from,]==which(clus_label_mat[from,]==q) & membership_mat[to,]==which(clus_label_mat[to,]==k))
    }
  }
  # N_nodes_mat = N_nodes_mat/rowSums(N_nodes_mat)
  N_nodes_mat_collapse = N_nodes_mat_collapse + N_nodes_mat
  fields::image.plot(N_nodes_mat,zlim=c(0,1), main=paste0(rho_vec[from]," to ",rho_vec[to]), xaxt="n",yaxt="n")
  text(x = seq(0,1,length.out = 4), y = -0.25, labels = 1:4, xpd = TRUE)
  text(y = seq(0,1,length.out = 4), x = -0.25, labels = 1:4, xpd = TRUE)
}

fields::image.plot(N_nodes_mat_collapse/rowSums(N_nodes_mat_collapse),zlim=c(0,1), main=paste0(""), xaxt="n",yaxt="n")
text(x = seq(0,1,length.out = 4), y = -0.25, labels = 1:4, xpd = TRUE)
text(y = seq(0,1,length.out = 4), x = -0.25, labels = 1:4, xpd = TRUE)
text(x = 0.5, y = -0.35, labels = "From", xpd = TRUE)
text(y = 0.5, x = -0.35, labels = "To", xpd = TRUE)


# SBM schematic -----------------------------------------------------------


col_pal = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral"))
set.seed(831)
prob_mat = matrix(c(0.9,0.2,0.2,0.6),nrow=2)
prob_mat[1,2] = prob_mat[2,1]
image(prob_mat,zlim=c(0,3),col=col_pal,xaxt="n",yaxt="n")

n_node = 30; n1=20; n2 = n_node-n1
mem_mat = cbind(c(rep(1,n1),rep(0,n2)),c(rep(0,n1),rep(1,n2)))

adj_mat = matrix(rbinom(n_node^2, 1, mem_mat %*% prob_mat %*% t(mem_mat)), nrow=n_node,ncol=n_node)
adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
image(adj_mat, zlim=c(0,3),col=col_pal,xaxt="n",yaxt="n",asp=1,bty="n", box=FALSE)


adj_mat = matrix(rnorm(n_node^2, mean=mem_mat %*% prob_mat %*% t(mem_mat), sd=0.1), nrow=n_node,ncol=n_node)
adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
image(adj_mat, zlim=c(0,3),col=col_pal,xaxt="n",yaxt="n",asp=1,bty="n", box=FALSE)

v_vec = runif(n_node, 0,2)
v_mat = matrix(nrow=n_node,ncol=n_node)
for (i in 1:n_node) {
  for (j in 1:n_node) {
    v_mat[i,j] = max(v_vec[i],v_vec[j])
  }
}
image(adj_mat+v_mat,zlim=c(0,3), col=col_pal,xaxt="n",yaxt="n",asp=1,bty="n", box=FALSE)
image(v_mat,zlim=c(0,3), col=col_pal,xaxt="n",yaxt="n",asp=1,bty="n", box=FALSE)


image(t(mem_mat),col=gray.colors(2,start = 1,end = 0.3),xaxt="n",yaxt="n")






# explain non-identifiability issue ---------------------------------------

t_vec = seq(0,100,0.05)
y_lim = c(0,0.2)
c1 = 40; c2 = 10
pdf_true_array = array(dim=c(2,2,length(t_vec)))
pdf_true_array[1,1,] = dnorm(t_vec,60-c1,3)
pdf_true_array[1,2,] = pdf_true_array[2,1,] = dnorm(t_vec,70-c1,5)
pdf_true_array[2,2,] = dgamma(t_vec,((60-c2)^2)/(10^2),(60-c2)/(10^2))


pdf_array_list = pdf_true_array
pdf_array_list[1,1,] = dnorm(t_vec,60,3)
pdf_array_list[1,2,] = pdf_array_list[2,1,] = dnorm(t_vec,70,5)
pdf_array_list[2,2,] = dgamma(t_vec,((60)^2)/(10^2),(60)/(10^2))


plot_pdf_array(pdf_array_list = pdf_array_list,pdf_true_array = pdf_true_array,t_vec = t_vec,y_lim = y_lim)


# Over-clustering -------------------------------------------------

r = results2[[29]]
edge_time_mat = r$network$edge_time_mat
t_vec = r$network$t_vec
bw = 1
membership_true = r$network$membership_true


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




get_pairwise_dist = function(edge_time_mat, clusters, t_vec=seq(0, 50, 0.05), bw=1)
{
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                       n0_vec = n0_vec, t_vec = t_vec, bw = bw)
  # need modification
  # degree_mat = get_node_degree_mat(edge_time_mat = edge_time_mat, clusters = clusters)
  dist_mat = pairwise_dist_mat(pdf_array = node_pdf_array, degree_mat = degree_mat, t_unit = 0.05)$dist_mat
  return(dist_mat)
}

dist_mat_iter_1 =  get_pairwise_dist(edge_time_mat = edge_time_mat, clusters = clusters_iter_1)
plot_tsne(dist_mat_iter_1, membership_true, membership_iter_2)
plot_tsne(dist_mat_iter_1, membership_true, membership_merged_2)



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






