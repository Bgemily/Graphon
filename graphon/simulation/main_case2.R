
# Case 2 ------------------------------------------------------------------

rm(list=ls())
source('./cluster_point_proc.R')

NSim = 1000
Ncores = 10

step_size = 0.02
SEED_vec = seq(189,1107,length.out=NSim)
results2 = vector("list", 0)
count = 1

library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results2 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  main(case=2, SEED, k=3, step_size = step_size)
}

# for (SEED in SEED_vec) {
#   print('=============')
#   cat('case2, trial:',count, '\n')
#   print('=============')
#   count = count+1
#   results2[[as.character(SEED)]]=main2(SEED, k=3, step_size = step_size)
# }

now = format(Sys.time(), "%Y%m%d_%H%M")
save.image(paste0('case2_NSim', NSim, '_', now, '.Rdata'))



