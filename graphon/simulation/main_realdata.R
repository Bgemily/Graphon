#!/usr/bin/env Rscript

# Real data ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


library("optparse")

option_list = list(
  make_option(c('-f', '--data'), type="character", default="EdgeTime_20150410.csv"),
  make_option(c("-k", "--N_clus"), type="integer", default=3),
  make_option(c("-m", "--MaxIter"), type="integer", default=10),
  make_option(c('-b', "--bw"), type="integer", default=10)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data = opt$data
data_path = paste0("./realdata/", data)
N_clus = opt$N_clus
MaxIter = opt$MaxIter
bw = opt$bw


Ncores = 20

# step_size = 0.02


library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

library(data.table)
data_folder = "./processed_FunctionalData/"
path.list=list.files(data_folder);

for(k in 1:length(path.list)){
  path=path.list[[k]]
  edge_time_mat=as.matrix(fread(paste(data_folder, path, '/EdgeTime.csv', sep='')))
  edge_time_mat = edge_time_mat[,-1]
  
  max_time = max(edge_time_mat[which(edge_time_mat<Inf)])
  t_vec = seq(0, max_time, length.out = 1000)
  
  result = list()
  result$clus_result = do_cluster(edge_time_mat = edge_time_mat, N_clus = N_clus, t_vec = t_vec, MaxIter = MaxIter, bw = bw)
  result$network = list(t_vec = t_vec, edge_time_mat = edge_time_mat)
  
  now = format(Sys.time(), "%Y%m%d_%H%M")
  save.image(paste0(path, '_Nclus', N_clus, '_MaxIter', MaxIter, '_', now, '.Rdata'))
  
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

