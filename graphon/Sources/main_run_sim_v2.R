#!/usr/bin/env Rscript


# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

 

# Simulation setup --------------------------------------------------------


library("optparse")

option_list = list(
  make_option(c("-n", "--N_trial"), type="integer", default=50, 
              help="number of repeated trials"),
  make_option("--split", type="integer", default=5)
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

N_trial_total = opt$N_trial
split = opt$split

N_trial = N_trial_total/split



### Parallel computing setup
library(foreach)
library(doParallel)
N_cores = 20
registerDoParallel(cores=N_cores)


### Common parameters
N_clus = 3
time_shift_struc = max
total_time = 200
const = 40
max_iter = 3 ### number of iters when updating time shifts


### Simulation-specific parameters
N_subj = 5 
N_subj_list = as.list(c(1,5,10))
N_node_vec = rep(90, N_subj)

clus_size_mat = matrix(30, nrow=5, ncol=N_clus)
time_shift_mean_vec = rep(10,N_clus) 
time_shift_rad = min(time_shift_mean_vec)
conn_patt_var = 1 ### alpha
conn_patt_sep = 1.3 ### beta
conn_prob_mean = 0.7  
conn_prob_rad = 0 


clus_size_mat_list = list(matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus))
time_shift_mean_vec_list = list(rep(10,N_clus), rep(20,N_clus), rep(30,N_clus), rep(40,N_clus), rep(50,N_clus))
conn_patt_var_list = list(0.5,1,2,3,4,6,8)
conn_patt_sep_list = list(1.6, 1.5, 1.4, 1.3, 1.2, 1.1)
conn_prob_mean_list = list(0.7,0.5,0.4,0.3,0.2,0.19,0.18,0.17,0.16,0.15,0.14,0.13)


# clus_size_vec_list = list(c(30,30,30), c(32,32,26), c(34,34,22), c(36,36,18), 
#                           c(38,38,14), c(40,40,10), c(42,42,6))






# Run simulations ---------------------------------------------------------

results1 = list()

now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(paste0("../Results/Rdata/", now_sim), recursive = TRUE, showWarnings = FALSE)
setwd(paste0("../Results/Rdata/", now_sim))
for (. in 1:split) {

  for (i in 1:length(N_subj_list)) {
    N_subj = N_subj_list[[i]]
    tmp <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e5,1)
      tryCatch(main_v2(SEED=SEED, N_subj=N_subj, conn_prob_mean=conn_prob_mean),
               error = function(x) print(SEED))
    }
    results1[[paste0("N_subj_",N_subj)]] = tmp
  }

  now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
  save.image(paste0('N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))

}

setwd(paste0("../../../Sources/"))




# Visualize results -------------------------------------------------------

file_list = list.files(path=paste0("../Results/Rdata/", now_sim), full.names = TRUE)

Lambda_mean_sq_err = extract_measurement(file_list = file_list, measurement = "Lambda_mean_sq_err")
ARI_mean = extract_measurement(file_list = file_list, measurement = "ARI_mean")
F_mean_sq_err = extract_measurement(file_list = file_list, measurement = "F_mean_sq_err")
v_mean_sq_err = extract_measurement(file_list = file_list, measurement = "v_mean_sq_err")
F_shape_mean_sq_err = extract_measurement(file_list = file_list, measurement = "F_shape_mean_sq_err")


dir.create(paste0("../Results/Plots/", now_sim), recursive = TRUE, showWarnings = FALSE)
setwd(paste0("../Results/Plots/", now_sim))

pdf(file = "Lambda_mse.pdf", width=4, height=4)
plot_jitter_boxplot(data = Lambda_mean_sq_err, ylab="Mean squared error of lambda")
dev.off()


pdf(file = "F_mse.pdf", width=4, height=4)
plot_jitter_boxplot(data = F_mean_sq_err, ylab="Mean squared error of F")
dev.off()

pdf(file = "F_shape_mse.pdf", width=4, height=4)
plot_jitter_boxplot(data = F_shape_mean_sq_err, ylab="Mean squared error of normalized F")
dev.off()


pdf(file = "Z_ARI.pdf", width=4, height=4)
plot_jitter_boxplot(data = ARI_mean, ylim=c(0,1), ylab="Mean ARI")
dev.off()

pdf(file = "V_mse.pdf", width=4, height=4)
plot_jitter_boxplot(data = v_mean_sq_err, ylab="Mean squared error of v")
dev.off()

setwd(paste0("../../../Sources/"))



# Test --------------------------------------------------------------------

# library(foreach)
# library(doParallel)
# N_cores = 1
# registerDoParallel(cores=N_cores)
# 
# split = 1
# N_trial = 1
# 
# results1 = list()
# 
# now = format(Sys.time(), "%Y%m%d_%H%M%S")
# dir.create(paste0("../Results/Rdata/", now), recursive = TRUE, showWarnings = FALSE)
# setwd(paste0("../Results/Rdata/", now))
# for (. in 1:split) {
#   
#   for (i in 1:length(N_subj_list)) {
#     # tmp <- foreach(j = 1:N_trial) %dopar% {
#     #   SEED = sample(1:1e5,1)
#     #   N_subj = N_subj_list[[i]]
#     #   tryCatch(main_v2(SEED=SEED, N_subj=N_subj),
#     #            error = function(x) print(SEED))
#     # }
#     SEED = sample(1:1e5,1)
#     N_subj = N_subj_list[[i]]
#     tmp = tryCatch(main_v2(SEED=SEED, N_subj=N_subj),
#              error = function(x) print(SEED))
#     results1[[paste0("N_subj_",N_subj)]] = tmp
#   }
#   
#   now = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0('N_trial', N_trial, '_our', '_', now, '.Rdata'))
#   
# }
# 
# setwd(paste0("../../../Sources/"))


