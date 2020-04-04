#!/usr/bin/env Rscript

# Case 2 ------------------------------------------------------------------

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
  make_option(c("-s", "--specc"), type="logical", default=FALSE, 
              help="TRUE for spectral clustering and FALSE for kmeans [default=%default]")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


NSim = opt$NSim
pp = opt$pp
specc = opt$specc

# pp = FALSE
# NSim = 1000
Ncores = 20

step_size = 0.02
SEED_vec = seq(189,1107,length.out=NSim)
results2 = vector("list", 0)
count = 1

library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results2 <- foreach(i = 1:NSim) %dopar%
{
  SEED = SEED_vec[i]
  if (specc) main_specc(case=2, SEED, k=3, step_size = step_size, pp=TRUE, h=1)
  else if (pp) main_pp(case=2, SEED, k=3, step_size = step_size, h=1)
  else main(case=2, SEED, k=3, step_size = step_size)
}

results2 = combn_subj_recluster(results2, group_size = 1, MaxIter = 2)

# main(case=2, SEED, k=3, step_size = step_size)->r
# main_pp(case=2, SEED, k=3, step_size = step_size, h=1)->r

# for (SEED in SEED_vec) {
#   print('=============')
#   cat('case2, trial:',count, '\n')
#   print('=============')
#   count = count+1
#   results2[[as.character(SEED)]]=main2(SEED, k=3, step_size = step_size)
# }

if(pp) {results2_pp = results2; rm(results2)}

now = format(Sys.time(), "%Y%m%d_%H%M")
{
  if (specc) save.image(paste0('specc_case2_NSim', NSim, '_', now, '.Rdata'))
  else if (!pp) save.image(paste0('case2_NSim', NSim, '_', now, '.Rdata'))
  else save.image(paste0('pp_case2_NSim', NSim, '_', now, '.Rdata'))
}