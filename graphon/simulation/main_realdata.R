#!/usr/bin/env Rscript

# Real data ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


library("optparse")

option_list = list(
  make_option(c('-f', '--data_folder'), type="character", default="./processed_FunctionalData/"),
  make_option(c("-k", "--N_clus"), type="integer", default=3),
  make_option(c("-m", "--MaxIter"), type="integer", default=10),
  make_option(c('-b', "--bw"), type="integer", default=5)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_folder = opt$data_folder
N_clus = opt$N_clus
MaxIter = opt$MaxIter
bw = opt$bw
standardize = F
step_size=0.04


Ncores = 20

# step_size = 0.02


library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

library(data.table)
# data_folder = "../processed_FunctionalData/"

window_length_vec = seq(480,40,-40) # =seq(10s,120s,10s)
rho_vec = seq(0.8,0.3,-0.1)

for (window_length in window_length_vec) {
  for (rho in rho_vec) {
    
    path.list=list.files(data_folder);
    # for(k in 1:length(path.list)){
    for(k in 1:2){
      path=path.list[[k]]
      
      # estimate and save edge_time_mat
      res = est_edge_time(window_length = window_length, window_step = window_length, 
                          rho = rho, processed_data_folder = data_folder, path = path)
      cor.full.ave = res$cor.full.ave
      # load(paste0(data_folder,path,'/cor_full_ave','_win',window_length,'_rho',rho,'.rdata'))
      # edge_time_mat = res$edge_time_mat
      
      # apply our method and save clustering result
      edge_time_mat=as.matrix(read.csv(paste(data_folder, path, '/EdgeTime.csv', sep='')))
      edge_time_mat = edge_time_mat[,-1]
      
      avai.inds = as.matrix(read.csv(paste(data_folder,path,'/AvaiNeurons.csv',sep='')))
      avai.inds=avai.inds[,-1];
      
      locs.all = as.matrix(read.csv(paste(data_folder,path,'/locs_all.csv',sep='')))
      locs.all = locs.all[,-1]
      locs=locs.all[avai.inds,]
      
      mnx.all = as.matrix(read.csv(paste(data_folder,path,'/mnx.csv',sep='')))
      mnx.all = mnx.all[,-1]
      mnx = mnx.all[avai.inds]
      
      
      max_time = max(edge_time_mat[which(edge_time_mat<Inf)])
      t_vec = seq(0, max_time+10, length.out = 1000)
      
      L_side = which(locs[,2]<0 & rowSums(edge_time_mat[,which(locs[,2]<0)]<Inf)>=2)
      R_side = which(locs[,2]>0 & rowSums(edge_time_mat[,which(locs[,2]>0)]<Inf)>=2)
      
      L_result = list()
      L_result$clus_result = do_cluster(edge_time_mat = edge_time_mat[L_side,L_side], N_clus = N_clus, 
                                        t_vec = t_vec, MaxIter = MaxIter, bw = bw, 
                                        standardize = standardize, step_size = step_size)
      L_result$network = list(t_vec = t_vec, edge_time_mat = edge_time_mat[L_side,L_side])
      L_result$id=L_side
      
      R_result = list()
      R_result$clus_result = do_cluster(edge_time_mat = edge_time_mat[R_side,R_side], N_clus = N_clus, 
                                        t_vec = t_vec, MaxIter = MaxIter, bw = bw, 
                                        standardize = standardize, step_size = step_size)
      R_result$network = list(t_vec = t_vec, edge_time_mat = edge_time_mat[R_side,R_side])
      R_result$id=R_side
      
      
      membership = numeric(nrow(locs))
      membership[L_side] = clus2mem(L_result$clus_result$clusters)
      membership[R_side] = clus2mem(R_result$clus_result$clusters)
      
      
      now = format(Sys.time(), "%Y%m%d_%H%M")
      save(list=c("edge_time_mat",'locs','mnx','L_result','R_result','membership'),file=paste0(path, '_Nclus', N_clus, "_win", window_length, "_rho",rho, '_now', now, '.Rdata'))
      
      
      # make plots
      visual_realdata(L_result, path, suffix=paste0('_Nclus', N_clus, "_win", window_length, 
                                                    "_rho",rho, '_now', now,"_L.pdf"))
      visual_realdata(R_result, path, suffix=paste0('_Nclus', N_clus, "_win", window_length, 
                                                    "_rho",rho, '_now', now,"_R.pdf"))
      
      dat.dFF=as.matrix(fread(paste(data_folder,path,'/dFF.csv',sep='')))
      dat.dFF=dat.dFF[,-1]
      na.inds= is.na(dat.dFF[,1]);
      reduced.dFF=dat.dFF[!na.inds,];
      visual_realdata_2(result = L_result, path = path, suffix=paste0('_Nclus', N_clus, "_win", window_length, 
                                                                      "_rho",rho, '_now', now,"_R.pdf"), 
                        edge.time = edge.time, locs = locs, reduced.dFF = reduced.dFF, 
                        window_length = window_length, window_step = window_length, 
                        cor.full.ave = cor.full.ave, member.ship = membership)
      visual_realdata_2(result = R_result, path = path, suffix=paste0('_Nclus', N_clus, "_win", window_length, 
                                                      "_rho",rho, '_now', now,"_R.pdf"), 
                        edge.time = edge.time, locs = locs, reduced.dFF = reduced.dFF, 
                        window_length = window_length, window_step = window_length, 
                        cor.full.ave = cor.full.ave, member.ship = membership)
      
      
      # write.csv(membership, paste('../processed_FunctionalData/',path,'/MembShip.csv',sep=''),col.names=F)
    }
    
  }
}


# edge_time_mat = as.matrix(read.csv(data_path))
# edge_time_mat = edge_time_mat[,-1]
# max_time = max(edge_time_mat[which(edge_time_mat<Inf)])
# t_vec = seq(0, max_time+20, length.out = 1000)
# 
# result = list()
# result$clus_result = do_cluster(edge_time_mat = edge_time_mat, N_clus = N_clus, t_vec = t_vec, MaxIter = MaxIter, bw = 5)
# result$network = list(t_vec = t_vec, edge_time_mat = edge_time_mat)
# 
# now = format(Sys.time(), "%Y%m%d_%H%M")
# save.image(paste0(data, '_Nclus', N_clus, '_MaxIter', MaxIter, '_', now, '.Rdata'))

