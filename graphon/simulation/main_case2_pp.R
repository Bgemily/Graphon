
# Case 2 ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

NSim = 100
Ncores = 10

step_size = 0.02
SEED_vec = seq(189,1107,length.out=NSim)
results2 = vector("list", 0)
count = 1
pp = TRUE

library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results2 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  if (pp) main_pp(case=2, SEED, k=3, step_size = step_size, h=1)
  else main(case=2, SEED, k=3, step_size = step_size)
}

# main(case=2, SEED, k=3, step_size = step_size)->r
# main_pp(case=2, SEED, k=3, step_size = step_size, h=1)->r

# for (SEED in SEED_vec) {
#   print('=============')
#   cat('case2, trial:',count, '\n')
#   print('=============')
#   count = count+1
#   results2[[as.character(SEED)]]=main2(SEED, k=3, step_size = step_size)
# }

now = format(Sys.time(), "%Y%m%d_%H%M")
{
  if (!pp) save.image(paste0('case2_NSim', NSim, '_', now, '.Rdata'))
  else save.image(paste0('pp_case2_NSim', NSim, '_', now, '.Rdata'))
}