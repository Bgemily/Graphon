# Case 3 ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

NSim = 100
Ncores = 10

step_size = 0.02
SEED_vec = seq(139,8397,length.out=NSim)
results3 = vector("list", 0)
count = 1
pp = TRUE

library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results3 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  if (pp) main_pp(case=3, SEED, k=3, step_size = step_size, h=1)
  else main(case=3, SEED, k=3, step_size = step_size)
}

# for (SEED in c(SEED_vec[881])) {
#   print('=============')
#   cat('case3, trial:',count, '\n')
#   print('=============')
#   count = count+1
#   results3[[as.character(SEED)]]=main3(SEED, k=3, step_size = step_size)
# }


now = format(Sys.time(), "%Y%m%d_%H%M")
{
  if (!pp) save.image(paste0('case3_NSim', NSim, '_', now, '.Rdata'))
  else save.image(paste0('pp_case3_NSim', NSim, '_', now, '.Rdata'))
}




