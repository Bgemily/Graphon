#!/usr/bin/env Rscript

# Case 1 ---------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)



library("optparse")

option_list = list(
  make_option(c("-n", "--NSim"), type="integer", default=1, 
              help="number of repeated trials")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


NSim = opt$NSim

Ncores = 10

# step_size = 0.02
SEED_vec = seq(1908,16417,length.out=NSim)


library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results1 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  main(case=1, SEED=SEED, N_clus=2, N_overclus=3, MaxIter = 2)
}


now = format(Sys.time(), "%Y%m%d_%H%M")
save.image(paste0('case1_NSim', NSim, '_', now, '.Rdata'))
  

