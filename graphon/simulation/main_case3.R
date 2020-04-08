#!/usr/bin/env Rscript

# Case 3 ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)



library("optparse")

option_list = list(
  make_option(c("-n", "--NSim"), type="integer", default=1, 
              help="number of repeated trials"),
  make_option(c("-p", "--pp"), type="logical", default=FALSE, 
              help="TRUE for pdf and FALSE for cdf [default=%default]"),
  make_option(c("-s", "--specc"), type="logical", default=TRUE, 
              help="TRUE for spectral clustering and FALSE for kmeans [default=%default]")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


NSim = opt$NSim
pp = opt$pp
specc = opt$specc


# pp = FALSE
# NSim = 1000
Ncores = 10

step_size = 0.02
SEED_vec = seq(139,8397,length.out=NSim)
results3 = vector("list", 0)
count = 1


library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results3 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  if (specc) main_specc(case=3, SEED, k=3, step_size = step_size, pp=TRUE, h=1)
  else if (pp) main_pp(case=3, SEED, k=3, step_size = step_size, h=1)
  else main(case=3, SEED, k=3, step_size = step_size)
}

results3 = combn_subj_recluster(results3, group_size = 1, MaxIter = 1)


# for (SEED in c(SEED_vec[881])) {
#   print('=============')
#   cat('case3, trial:',count, '\n')
#   print('=============')
#   count = count+1
#   results3[[as.character(SEED)]]=main3(SEED, k=3, step_size = step_size)
# }

if(pp) {results3_pp = results3; rm(results3)}

now = format(Sys.time(), "%Y%m%d_%H%M")
{
  if (specc) save.image(paste0('specc_case3_NSim', NSim, '_', now, '.Rdata'))
  else if (!pp) save.image(paste0('case3_NSim', NSim, '_', now, '.Rdata'))
  else save.image(paste0('pp_case3_NSim', NSim, '_', now, '.Rdata'))
}




