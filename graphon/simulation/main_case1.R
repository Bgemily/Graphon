# source('~/Documents/Academic/SC/graphon/simulation/cluster_point_proc.R')
# source('./cluster_point_proc.R')

# Case 1 ---------------------------------------------------------------

rm(list=ls())
source('./cluster_point_proc.R')

NSim = 1000
Ncores = 10

step_size = 0.02
SEED_vec = seq(1908,16417,length.out=NSim)
results = vector("list", 0)
count = 1

library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results1 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  main(case=1, SEED, k=2, step_size = step_size)
}


# for (SEED in SEED_vec) {
#   print('================')
#   cat('case1, trial:', count, '\n')
#   print('================')
#   count = count+1
#   results[[as.character(SEED)]]=main(SEED, k=2, step_size = step_size)
# }


now = format(Sys.time(), "%Y%m%d_%H%M")
save.image(paste0('case1_NSim', NSim, '_', now, '.Rdata'))



