

generate_network2_v2 = function(N_subj, N_node_vec, N_clus=3, total_time=200, t_vec=seq(0,total_time,length.out=1000), 
                                clus_size_mat=matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
                                conn_patt_var=1, conn_patt_sep=2, conn_prob_mean=0.8, conn_prob_rad=0.1, 
                                time_shift_mean_vec=rep(10,N_clus), time_shift_rad=min(time_shift_mean_vec), 
                                time_shift_struc=max, 
                                SEED=NULL, const=40)
{
  if(!is.null(SEED)) set.seed(SEED)
  
  if (sum(rowSums(clus_size_mat)!=N_node_vec)>0) {
    stop("N_node_vec is inconsistent with clus_size_mat.")
  }

# Generate cluster memberships ---------------------------------------

  membership_true_list = clus_true_list = vector(mode="list", length=N_subj)
  for (i in 1:N_subj) {
    membership_tmp = rep(1:N_clus, clus_size_mat[i,])
    membership_true_list[[i]] = membership_tmp
    clus_true_list[[i]] = mem2clus(membership_tmp)
  }
  

# Generate time shifts ----------------------------------------------------

  time_shift_list = vector(mode = "list", length = N_subj)
  
  if (length(time_shift_mean_vec) != N_clus) {
    stop("Length of time_shift_mean_vec should be equal to N_clus.")
  }
  
  time_shift_max_vec = time_shift_mean_vec + time_shift_rad
  time_shift_min_vec = time_shift_mean_vec - time_shift_rad
  
  for (m in 1:N_subj) {
    time_shift_tmp = sapply(1:N_clus, function(k)runif(n = clus_size_mat[m,k], 
                                                       min = time_shift_min_vec[k],
                                                       max = time_shift_max_vec[k]) )
    if (is.list(time_shift_tmp)) {
      time_shift_tmp = unlist(time_shift_tmp)
    }
    time_shift_tmp = c(time_shift_tmp)
    
    time_shift_list[[m]] = time_shift_tmp
  }
  


# Generate connecting pattern functions --------------------------------------------

  mean_mat = matrix(nrow = N_clus, ncol = N_clus)
  var_mat = matrix(nrow = N_clus, ncol = N_clus)

  mean_mat[1,2] = mean_mat[2,1] = const*conn_patt_sep
  var_mat[1,2] = var_mat[2,1] = const*conn_patt_var*conn_patt_sep^2
  
  mean_mat[1,3] = mean_mat[3,1] = const/conn_patt_sep^2
  var_mat[1,3] = var_mat[3,1] = const*conn_patt_var/conn_patt_sep
  
  mean_mat[2,3] = mean_mat[3,2] = const/conn_patt_sep 
  var_mat[2,3] = var_mat[3,2] = const*conn_patt_var/conn_patt_sep^2 
  
  mean_mat[1,1] = const
  var_mat[1,1] = const*conn_patt_var
  
  mean_mat[2,2] = const*conn_patt_sep^2
  var_mat[2,2] = const*conn_patt_var*conn_patt_sep
  
  mean_mat[3,3] = const*sqrt(conn_patt_sep)
  var_mat[3,3] = const*conn_patt_var*conn_patt_sep

  conn_prob_min = conn_prob_mean - conn_prob_rad
  conn_prob_max = conn_prob_mean + conn_prob_rad
  if (conn_prob_max > 1 | conn_prob_min < 0) {
    stop("Invalid conn_prob_rad.")
  }
  
  conn_prob_mat = matrix(ncol=N_clus, nrow=N_clus)
  conn_prob_mat[!lower.tri(conn_prob_mat)] = seq( from=conn_prob_min, to=conn_prob_max, length.out=((N_clus^2+N_clus)/2) )
  conn_prob_mat[lower.tri(conn_prob_mat)] = t(conn_prob_mat)[lower.tri(conn_prob_mat)]
  
  if (!isSymmetric(conn_prob_mat))
    stop("Constructed conn_prob_mat is not symmetric.")
  
  ### Specify connecting patterns' corresponding pdf functions and data generating functions
  # {  
  # pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  # pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=mean_mat[1,2]^2/var_mat[1,2], 
  #                                                               rate=mean_mat[1,2]/var_mat[1,2])*conn_prob_mat[1,2], 
  #                                      random = function(n) ifelse(test = runif(n)<=conn_prob_mat[1,2], 
  #                                                                  yes = rgamma(n, shape=mean_mat[1,2]^2/var_mat[1,2], 
  #                                                                               rate=mean_mat[1,2]/var_mat[1,2]), 
  #                                                                  no = Inf)); 
  # pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  # pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=mean_mat[1,3]^2/var_mat[1,3], 
  #                                                               rate=mean_mat[1,3]/var_mat[1,3])*conn_prob_mat[1,3], 
  #                                      random = function(n) ifelse(test = runif(n)<=conn_prob_mat[1,3], 
  #                                                                  yes = rgamma(n, shape=mean_mat[1,3]^2/var_mat[1,3], 
  #                                                                               rate=mean_mat[1,3]/var_mat[1,3]), 
  #                                                                  no = Inf));
  # pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]]
  # pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dgamma(x, shape=mean_mat[2,3]^2/var_mat[2,3], 
  #                                                               rate=mean_mat[2,3]/var_mat[2,3])*conn_prob_mat[2,3], 
  #                                      random = function(n) ifelse(test = runif(n)<=conn_prob_mat[2,3], 
  #                                                                  yes = rgamma(n, shape=mean_mat[2,3]^2/var_mat[2,3], 
  #                                                                               rate=mean_mat[2,3]/var_mat[2,3]), 
  #                                                                  no = Inf));
  # pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  # 
  # pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dgamma(x, shape=mean_mat[1,1]^2/var_mat[1,1], 
  #                                                               rate=mean_mat[1,1]/var_mat[1,1])*conn_prob_mat[1,1], 
  #                                      random = function(n) ifelse(test = runif(n)<=conn_prob_mat[1,1], 
  #                                                                  yes = rgamma(n, shape=mean_mat[1,1]^2/var_mat[1,1], 
  #                                                                               rate=mean_mat[1,1]/var_mat[1,1]), 
  #                                                                  no = Inf))
  # pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dgamma(x, shape=mean_mat[2,2]^2/var_mat[2,2], 
  #                                                               rate=mean_mat[2,2]/var_mat[2,2])*conn_prob_mat[2,2], 
  #                                      random = function(n) ifelse(test = runif(n)<=conn_prob_mat[2,2], 
  #                                                                  yes = rgamma(n, shape=mean_mat[2,2]^2/var_mat[2,2], 
  #                                                                         rate=mean_mat[2,2]/var_mat[2,2]), 
  #                                                                  no = Inf) )
  # pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dgamma(x, shape=mean_mat[3,3]^2/var_mat[3,3], 
  #                                                               rate=mean_mat[3,3]/var_mat[3,3])*conn_prob_mat[3,3], 
  #                                      random = function(n) ifelse(test = runif(n)<=conn_prob_mat[3,3], 
  #                                                                  yes = rgamma(n, shape=mean_mat[3,3]^2/var_mat[3,3], 
  #                                                                         rate=mean_mat[3,3]/var_mat[3,3]), 
  #                                                                  no = Inf))
  # 
  #   }
  

  ### Evaluate true connecting patterns
  pdf_true_array = cdf_true_array = array(dim = c(N_clus, N_clus, length(t_vec)))
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      pdf_true_array[q,l,] = dgamma(x=t_vec, shape=mean_mat[q,l]^2/var_mat[q,l], 
                                    rate=mean_mat[q,l]/var_mat[q,l]) * conn_prob_mat[q,l]
      
      cdf_true_array[q,l,] = pgamma(q=t_vec, shape=mean_mat[q,l]^2/var_mat[q,l], 
                                    rate=mean_mat[q,l]/var_mat[q,l]) * conn_prob_mat[q,l]
    }
  }
  
  

# Generate edge time matrices ---------------------------------------------

  edge_time_mat_list = list()
  for (m in 1:N_subj) {
    N_node_tmp = N_node_vec[[m]]
    
    edge_time_mat_tmp = matrix(Inf, nrow = N_node_tmp, ncol = N_node_tmp)
    for (q in 1:N_clus) {
      for (l in 1:N_clus) {
        samples = ifelse(test = runif(clus_size_mat[m,q]*clus_size_mat[m,l])<=conn_prob_mat[q,l], 
                         yes = rgamma(n = clus_size_mat[m,q]*clus_size_mat[m,l], 
                                      shape=mean_mat[q,l]^2/var_mat[q,l], 
                                      rate=mean_mat[q,l]/var_mat[q,l]), 
                         no = Inf)
        samples = matrix(samples, clus_size_mat[m,q], clus_size_mat[m,l]) 
        edge_time_mat_tmp[clus_true_list[[m]][[q]], clus_true_list[[m]][[l]]] = samples
      }
    }
    edge_time_mat_tmp[lower.tri(edge_time_mat_tmp)] = t(edge_time_mat_tmp)[lower.tri(edge_time_mat_tmp)] # make it symmetric
    
    ### Compute time shift matrix
    time_shift_vec_tmp = time_shift_list[[m]]
    time_shift_mat_tmp = matrix(nrow = N_node_tmp, ncol = N_node_tmp)
    if (is.function(time_shift_struc)){
      for (i in 1:N_node_tmp) {
        for (j in 1:N_node_tmp) {
          time_shift_mat_tmp[i,j] = time_shift_struc(time_shift_vec_tmp[i], time_shift_vec_tmp[j])
        }
      }
    }
    else
      stop("Invalid time_shift_struc. Supposed to be a function.")
    
    
    edge_time_mat_tmp = edge_time_mat_tmp + time_shift_mat_tmp # add time shifts 
    # edge_time_mat_tmp[ edge_time_mat_tmp>total_time ] = Inf
    
    edge_time_mat_list[[m]] = edge_time_mat_tmp
  }
  
  


# Output ------------------------------------------------------------------

  
  return(list(edge_time_mat_list=edge_time_mat_list, 
              membership_true_list=membership_true_list, clus_true_list=clus_true_list,
              time_shift_list=time_shift_list, 
              cdf_true_array=cdf_true_array, pdf_true_array=pdf_true_array
              ))
  
}


# Test --------------------------------------------------------------------


# res = generate_network2_v2(N_subj = 2,N_node_vec = c(30,30), N_clus = 3,total_time = 200,conn_patt_var = 2,
#                            conn_patt_sep = 2, conn_prob_mean = 0.8,conn_prob_rad = 0, 
#                            time_shift_mean_vec = rep(10,3),time_shift_rad = 5,)
# 
# 
# par(mfrow=c(1,1))
# image(res$edge_time_mat_list[[1]])
# par(mfrow=c(3,3))
# for (q in 1:3) {
#   for (k in 1:3) {
#     plot(res$pdf_true_array[q,k,])
#   }
# }
# res$time_shift_list
